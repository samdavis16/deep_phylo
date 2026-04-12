import atexit
import copy
import csv
import gzip
import os
import shutil
import socket
import sqlite3
import threading
import urllib.error
import urllib.request
from concurrent import futures
from pathlib import Path

import ete3
import matplotlib as mpl
import requests
import xmltodict
from Bio import Seq, SeqIO, SeqRecord
from ete3.parser.newick import NewickError

from . import file_util
from . import tree
from . import cluster


def clade_assignment(tree_file, boundaries):
    """ Returns a dictionary mapping ID to clade assignment on the basis of provide boundaries of the form
     [clade_name, [in_nodes], out_nodes ] where in_nodes infer a LCA of the clade. """

    try:
        t = ete3.Tree(tree_file)
    except NewickError:
        t = ete3.Tree(tree_file, format=2)

    assign = {}
    for clade in boundaries:
        t.set_outgroup(t&clade[2])  # Root on out_node for given clade
        leaf_names = tree.get_subtree_leaves(t, ext_nodes=clade[1])
        for name in leaf_names:
            assign[name] = clade[0]

    return assign


def custom_filter(
        seqs,
        function,
        pass_value=True,
        out_file=None,
        no_return=False
):
    """ Filter out sequences which fail a given test, i.e. do not return the required pass_value (default=True) for the
    provided lambda function. """

    if no_return and not out_file:
        raise RuntimeError("If no return, an output file name must be specified.")

    if isinstance(seqs, str):
        seqs = list(SeqIO.parse(seqs, 'fasta'))

    retain_seqs = [seq for seq in seqs if function(seq) == pass_value]

    if out_file:
        SeqIO.write(retain_seqs, out_file, 'fasta')

    if not no_return:
        return retain_seqs


def fetch_seq_record(acc, dbName='uniprotkb'):
    """ Return a Bio.SeqRecord object from a given accession and database.
     TODO: Currently only UniprotKB handled here. """

    raw_annots = fetch_annot(acc, dbName=dbName)

    if dbName == 'uniprotkb':

        if raw_annots:  # We got annot data
            seq = ""
            seq_idx = [line.startswith("SQ") for line in raw_annots].index(True)  # Where does the sequence info begin?
            for i in range(seq_idx+1, len(raw_annots)-2):
                line = raw_annots[i]
                subseq = line.replace(' ', '')
                seq += subseq
            seq_record = SeqRecord.SeqRecord(Seq.Seq(seq), id=acc, name=acc)
            return seq_record

        else:  # Couldn't get annots so warn user and return None
            print(f"Sequence fetch failed for {acc}")
            return

    elif dbName == 'uniparc':

        # Construct URL
        url = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?style=raw&Retrieve=Retrieve&db=uniparc&format=fasta&id=' + acc

        try:
            data = urllib.request.urlopen(url).read().decode("utf-8")
            if data.startswith("ERROR"):
                return None  # Assume provided ID is not uniparc, so return None
            # Collect sequence from FASTA format
            seq = ''
            for line in data.split('\n')[1:]:
                seq += line
            return SeqRecord.SeqRecord(Seq.Seq(seq), id=acc, name=acc)
        except urllib.error.HTTPError as ex:  # Assume non-query related issue (e.g. connection issue). Warn user but keep going
            raise RuntimeWarning(ex.read())


def fetch_seq_records(accs, dbName='uniprotkb'):
    """ Return a list of Bio.SeqRecord objects from a list of accession numbers and selected databse. """

    records = []
    for acc in accs:
        seq_record = fetch_seq_record(acc, dbName)
        if seq_record:
            records.append(seq_record)
    return records


def fetch_annot_raw():
    return


# TODO: All fetching functions should call fetch_batch_raw
def fetch_batch_raw(
        entryIds,
        dbName='uniprotkb',
        format='uniprot',
        report_failures=True
):
    """ Fetch a batch of up to 200 queries from a given dbfetch database. Results are in a single with new lines
    separating. If report_failures is True, results are checked for the presence of each provided ID, and a double
    ('results', [missing_ids]) is returned. Else, only the results string is returned. """

    # Base url for all EBI web services
    __ebiUrl__ = 'http://www.ebi.ac.uk/Tools/'

    if isinstance(entryIds, str):
        entryIds = [entryIds]
    elif not isinstance(entryIds, list):
        raise RuntimeError("entryIds must be a string (single query) or list of strings (multiple queries).")

    missing = []  # Keep track of IDs which fail for fetching
    while entryIds:

        attempt_n = 1  # Allow 10 attempts for each version of batch list
        acc_str = ','.join(entryIds)
        url = __ebiUrl__ + 'dbfetch/dbfetch?style=raw&Retrieve=Retrieve&db=' + dbName + '&format=' + format + '&id=' + acc_str

        while True:
            try:
                data = urllib.request.urlopen(url, timeout=60).read().decode("utf-8")  # Query dbfetch
                break

            except urllib.error.HTTPError as ex:  # Assume non-query related issue (e.g. connection)

                # TODO: For debugging - DELETE
                # with open("annot_workflow_debug.log", "a") as error_file:
                #     error_file.write(f"{ex.read()}\tAttempt {attempt_n}\n")

                # TODO: For now, terminate after 10 unsuccessful attempts on same batch - need to implement better handling here
                print(f"{ex.read()}\tAttempt {attempt_n}")
                if attempt_n == 10:  # Terminate search
                    print(f"Max attempts reached: terminating for batch: {acc_str}")
                    # # TODO: For debugging - DELETE
                    # with open("annot_workflow_debug.log", "a") as error_file:
                    #     error_file.write(f"Max attempts reached: terminating for batch: {acc_str}\n")
                    for acc in entryIds:  # Add remaining accessions to missing list
                        missing.append(acc)
                    entryIds = []
                    data = ""  # Empty string is returned
                    break
                else:
                    # with open("annot_workflow_debug.log", "a") as error_file:
                    #     error_file.write(f"Failure on attempt {attempt_n} for {acc_str}\n")
                    attempt_n += 1

            except urllib.error.URLError as ex:
                if attempt_n == 10:  # Terminate search

                    for acc in entryIds:  # Add remaining accessions to missing list
                        missing.append(acc)
                    entryIds = []
                    data = ""  # Empty string is returned
                    break
                else:
                    attempt_n += 1

            except socket.timeout:
                if attempt_n == 10:  # Terminate search

                    for acc in entryIds:  # Add remaining accessions to missing list
                        missing.append(acc)
                    entryIds = []
                    data = ""  # Empty string is returned
                    break
                else:
                    attempt_n += 1

        # Check for presence of annots for all IDs
        to_remove = None
        for acc in entryIds:
            if acc not in data:  # Case where first ID is missing and error is returned should be handled by this
                to_remove = acc
                break
        if to_remove:
            missing.append(to_remove)
            entryIds.remove(to_remove)
        else:  # Complete success on the current batch or entryIds is empty
            break

    # If 0 entries fetched successfully, data str will begin with 'ERROR'
    if data.startswith("ERROR"):
        data = ""
        # # TODO: For debugging - DELETE
        # with open("annot_workflow_debug.log", "a") as error_file:
        #     error_file.write(f"Server returned ERROR for {acc_str}\n")

    # Return results
    if report_failures:
        return data, missing
    else:
        return data


def fetch_batch_to_file(entryIds, result_file, lock, dbName='uniprotkb', format = 'uniprot', log_file=None):
    """ Write query results and (optionally) error log for a given batch to master output file(s) in a threadsafe
     manner. """

    # Fetch annots
    if log_file:
        # print(f"BATCH STARTED: {entryIds}")  #TODO: For testing, delete
        batch_result, missing_seqs = fetch_batch_raw(entryIds, dbName, format, report_failures=True)
        # Write missing seqs to log separated by new lines
        if missing_seqs:
            batch_log = '\n'.join(missing_seqs) + '\n'
        else:
            batch_log = ''
    else:
        batch_result = fetch_batch_raw(entryIds, dbName, format)

    with lock:  # Wait to acquire lock
        # Write results and log (if provided)
        with open(result_file, 'a') as res_file:
            res_file.write(batch_result)
        if log_file:
            with open(log_file, 'a') as file:
                file.write(batch_log)
            # TODO: For debugging - DELETE
        with open("test_debug.log", 'a') as log:
            log.write("Finished a thread\n")


def fetch_raw_to_file(
        entryIds,
        outFile,
        dbName='uniprotkb',
        format='uniprot',
        failure_log=None,
        max_batch_size=200,
        max_threads=10
):
    """ Fetch raw annotations for a given dbfetch database for a single query ID or list. If a failure log file is
    specified, queries not successfully retrieved are reported. Output file is all raw annotation strings concatenated.
    Batch size default is dbfetch default (200). For sets of IDs where a high rate of fetch failure is expected, use
    lower batch size and more threads. NOTE: Uniparc fetching only allows 25 per batch for some reason. """

    if isinstance(entryIds, str):
        entryIds = [entryIds]
    elif not isinstance(entryIds, list):
        raise RuntimeError("entryIds must be a string (single query) or list of strings (multiple queries).")

    # Make sure max batch size is not more than 25 for uniparc!
    if dbName == 'uniparc' and max_batch_size > 25:
        max_batch_size = 25

    # Divide into batches (default max is dbfetch batch max)
    batches = []
    for i in range(0, len(entryIds), max_batch_size):
        batches.append(entryIds[i:i+max_batch_size])

    # TODO: For debugging - DELETE
    with open("test_debug.log", 'a') as log:
        log.write(f"N batches: {len(batches)}\n")

    # # TODO: For testing, remove (batch validation)
    # already_seen = []
    # for batch in batches:
    #     already_this = []
    #     for acc in batch:
    #         if acc in already_this:
    #             print(f"Intra-batch duplicate: {acc}")
    #         else:
    #             already_this.append(acc)
    #             if acc in already_seen:
    #                 print(f"Multiple batches: {acc}")
    #             else:
    #                 already_seen.append(acc)
    # for acc in entryIds:
    #     if acc not in already_seen:
    #         print(f"Missing: {acc}")

    # Delete files of same name if they exist
    if os.path.exists(outFile):
        os.remove(outFile)

    if failure_log:
        if os.path.exists(failure_log):
            os.remove(failure_log)

    executor = futures.ThreadPoolExecutor(max_workers=max_threads)
    lock = threading.Lock()

    # Execute batch fetching by threads
    future_store = [executor.submit(fetch_batch_to_file, batches[i], outFile, lock, dbName, format, failure_log)
                    for i in range(len(batches))]

    # TODO: For debugging - DELETE
    with open("test_debug.log", 'a') as log:
        log.write("Completed threads\n")

    atexit.register(executor.shutdown, wait=False)

    futures.wait(future_store)


def fetch_uparc_taxonomy(
        uparc_annots,
        tax_file,
        tax_idx=None,
        uparc_idx=None,
        uparc_ids=None
):
    """ Fetch NCBI taxonomy information for all or selected entries in a Uniparc annotation file. """

    if not uparc_idx:  # Need to make index file
        uparc_idx = uparc_annots.split('.')[0] + ".idx"
        if uparc_idx not in os.listdir():
            uniparc_annot_idx(uparc_annots, uparc_idx)

    if not tax_idx:
        tax_idx = tax_file.split(".")[0] + ".idx"

    idx_map = file_util.tab_del_file_to_dict(uparc_idx)

    if not uparc_ids:  # Need to fetch for all entries
        to_fetch = {acc : idx for acc, idx in idx_map.items()}
    else:
        to_fetch = {acc : idx_map[acc] for acc in uparc_ids if acc in idx_map.keys()}

    tax_ids = {}  # Map Uniparc ID to taxonomy ID
    for acc, idx in to_fetch.items():
        xml_dict = xmltodict.parse(uparc_annots_from_file(uparc_annots, idx=int(idx)))
        descriptions = xml_dict["rdf:RDF"]["rdf:Description"]
        tax_id = None
        for description in descriptions:
            if "organism" in description.keys():
                tax_url = description["organism"]["@rdf:resource"]
                if len(tax_url.split("/")) > 1 and tax_url.split("/")[-2] == "taxonomy":
                    tax_id = tax_url.split("/")[-1]
                    break
        if tax_id:
            tax_ids[acc] = tax_id

    # Collect all NCBI taxonomy record for sequences for which an acc was mapped
    fetch_raw_to_file(list(tax_ids.values()), tax_file, dbName='taxonomy', format='enataxonomyxml',
                      failure_log=tax_file.split(".")[0]+".log", max_batch_size=200,
                      max_threads=10)

    # Taxonomy index file is compulsory, otherwise link to the Uniparc accession is lost
    with open(tax_file) as annot_file:
        with open(tax_idx if tax_idx else tax_file.split(".")[0]+".idx", 'w') as idx_file:
            pos = 0
            allocated = []  # Currently tax file may have duplicate entries, avoid matching to Uparc accs multiple times
            for line in annot_file:
                if line.startswith("<taxon") and "rank=" in line:  # Found the start of a tax record
                    this_tax_id = line.split('taxId="')[1].split('"')[0]
                    if this_tax_id not in allocated:
                        for acc, tax_id in tax_ids.items():
                            if this_tax_id == tax_id:
                                idx_file.write(f"{acc}\t{pos}\n")
                        allocated.append(this_tax_id)
                pos += len(line)


def uniprot_annot_idx(raw_annot_file, idx_file):
    """Create offset index for raw Uniprot annotation file."""
    pos = 0
    with open(raw_annot_file) as annots:
        with open(idx_file, 'w') as offset_file:
            while True:
                line = annots.readline()
                if not line:
                    break
                if line.startswith("ID"):
                    offset = pos
                    pos += len(line)
                    line = annots.readline()  # Read accession line
                    acc = line.strip().split()[1].split(';')[0]
                    offset_file.write(f"{acc}\t{offset}\n")
                pos += len(line)


def uniparc_annot_idx(raw_annot_file, idx_file):
    """ Create offset index for raw (xml format) Uniparc annotation file. Entries contain the <entry> object. Any
     outer objects which may represent batches of search results which have been subsequently concatenated are ignored.
     """

    with open("test_debug.log", 'a') as log:
        log.write("Started annot_idx\n")

    pos = 0
    with open(raw_annot_file) as annots:
        with open(idx_file, 'w') as offset_file:
            while True:
                line = annots.readline()
                if not line:  # EOF
                    break
                if line.startswith("<rdf:RDF"):  # Beginning of Uniparc entry record
                    offset = pos  # Beginning position
                elif line.startswith('<rdf:Description rdf:about="UPI') and "#" not in line and "faldo" not in line:  # Found Uniparc ID
                    acc = line.split('"')[1]
                    offset_file.write(f"{acc}\t{offset}\n")  # Now have acc and offset so write to index file
                pos += len(line)


def merge_raw_annot_files(
        in_files,
        out_file,
        idx_files="present",
        record_type='uniprot',
        remove_input_files=False
):
    """ Merge multiple raw annotation files into a single one, removing any duplicate records. If index files are
    already present and obeying naming conventions (.idx for the corresponding .txt annotation file), use default flag.
    If files are present but have different names, provide as a list. If indexes are not present, set idx_files=None. """

    if not isinstance(in_files, list):
        raise RuntimeError("At least two files must be provided for merging.")

    if out_file in in_files:  # If an existing filename is to be merged into
        out_name = f"temp_{out_file}"
    else:
        out_name = out_file

    if not idx_files:  # If index files aren't provided, generate them

        in_idxs = {}
        for file in in_files:
            idx_name = file.split('.txt')[0]+'.idx'
            if record_type == 'uniprot':
                uniprot_annot_idx(file, idx_name)
            elif record_type == 'uniparc':
                uniparc_annot_idx(file, idx_name)
            else:
                raise RuntimeError(f"{record_type} is not a supported annotation type.")
            in_idxs[file] = file_util.tab_del_file_to_dict(idx_name)

    elif idx_files == "present":   # Index files are present and obey naming conventions
        in_idxs = {file : file_util.tab_del_file_to_dict(file.split('.txt')[0] + '.idx') for file in in_files}

    elif isinstance(idx_files, list):   # Index files are present but named differently
        in_idxs = {in_files[i] : file_util.tab_del_file_to_dict(idx_files[i]) for i in range(len(idx_files))}

    else:
        raise RuntimeError("Invalid idx_file parameter provided.")

    # Get non-redundant set of sequence IDs
    seq_ids = set()
    for idx_dict in in_idxs.values():
        seq_ids.update(set(idx_dict.keys()))

    # Write annotation for each sequence once to merged file
    with open(out_name, 'w') as file:
        for seq in seq_ids:
            for in_file, idx in in_idxs.items():
                if seq in idx.keys():
                    if record_type == 'uniprot':
                        file.write(uprot_annots_from_file(in_file, idx=idx[seq]))
                    elif record_type == 'uniparc':
                        file.write(uparc_annots_from_file(in_file, idx=idx[seq]))
                    break

    # TODO: Uncomment below once confirmed that everything is working - also change next statement from if to elif
    # if remove_input_files:  # Remove old files if req
    #     for i in range(len(in_files)):
    #         os.remove(in_files[i])
    #         if isinstance(idx_files, list):
    #             os.remove(idx_files[i])
    #         else:
    #             os.remove(idx_files[i].split('.txt')[0]+'.idx')

    if out_name == f"temp_{out_file}":  # if merging all into existing file name, rename
        os.rename(out_name, out_file)

    # Create index for merged annotations
    if record_type == 'uniprot':
        uniprot_annot_idx(out_file, out_file.split('.txt')[0]+'.idx')
    elif record_type == 'uniparc':
        uniparc_annot_idx(out_file, out_file.split('.txt')[0]+'.idx')


def up_id_mapping(
        raw_uniparc,
        accs=None,
        parc2prot_file=None,
        prot2parc_file=None,
        idx_file=None,
        report_missing=True
):
    """
         Create mapping of Uniparc to UniprotKB (and/or vice-versa) IDs from a raw Uniparc annotation file (XML). Output
        file(s) are tab-delimited. For Uniparc IDs with multiple mapping UniprotKB IDs, IDs are listed separated with
        commas. If specific accessions are not provided, all annotations in the file are mapped.

        #TODO: Currently looking for missing IDs in index file but not annotation file """

    if not (parc2prot_file or prot2parc_file):
        raise RuntimeError("At least one input file must be specified for ID mapping.")

    if not idx_file:

        idx_file = raw_uniparc.split(".")[-1]+".idx"

        # Check if an idx file with specified naming format is available in working directory
        if idx_file not in os.listdir():
            uniparc_annot_idx(raw_uniparc, idx_file)

    mappings = {}  # Store mappings as { uparc_acc : [uprot_accs] }

    idx_map = file_util.tab_del_file_to_dict(idx_file)

    if accs:  # Only map specific requested uparc IDs
        idxs = []
        missing_ids = []
        for acc in accs:
            try:
                idxs.append((acc, int(idx_map[acc])))
            except KeyError:
                missing_ids.append(acc)

    else:  # We want all
        idxs = [(acc, int(idx)) for acc, idx in idx_map.items()]

    idxs.sort(key=lambda x: x[1])

    for acc, idx in idxs:  # Iterate through

        # Extract xml text for this entry, load in as xml dict and get dictionary of mapped source database entries
        db_sources = extract_uparc_sources(uparc_annots_from_file(raw_uniparc, idx=idx))

        if db_sources and "uniprot" in db_sources.keys():
            mappings[acc] = db_sources["uniprot"]
        else:
            mappings[acc] = []

    if accs and report_missing:
        print(f"The following IDs were not found in the provided index file: {', '.join(missing_ids)}")

    # Write uparc2uprot file if specified
    if parc2prot_file:
        with open(parc2prot_file, 'w') as out_file:
            for uparc, uprots in mappings.items():
                out_file.write(f"{uparc}\t{','.join(uprots)}\n")

    # Write uprot2uparc file if specified
    if prot2parc_file:
        reverse_mappings = {}
        for uparc, uprots in mappings.items():
            for uprot in uprots:
                reverse_mappings[uprot] = uparc
        with open(prot2parc_file, 'w') as out_file:
            for uprot, uparc in reverse_mappings.items():
                out_file.write(f"{uprot}\t{uparc}\n")


def up_id_mapping_old(
        raw_uniparc,
        accs=None,
        parc2prot_file=None,
        prot2parc_file=None,
        idx_file=None,
        report_missing=True
):
    """
     NOTE: Function is no longer valid due to outdated Uniparc xml format. Updated function in testing.

     Create mapping of Uniparc to UniprotKB (and/or vice-versa) IDs from a raw Uniparc annotation file (XML). Output
    file(s) are tab-delimited. For Uniparc IDs with multiple mapping UniprotKB IDs, IDs are listed separated with
    commas. If specific accessions are not provided, all annotations in the file are mapped.

    #TODO: Currently looking for missing IDs in index file but not annotation file """

    if not (parc2prot_file or prot2parc_file):
        raise RuntimeError("At least one input file must be specified for ID mapping.")

    mappings = {}  # Store mappings as { uparc_acc : [uprot_accs] }

    # If index file available
    if idx_file:
        idx_map = file_util.tab_del_file_to_dict(idx_file)

        if accs:  # Only map specific requested uparc IDs
            idxs = []
            missing_ids = []
            for acc in accs:
                try:
                    idxs.append((acc, int(idx_map[acc])))
                except KeyError:
                    missing_ids.append(acc)

        else:  # We want all
            idxs = [(acc, int(idx)) for acc, idx in idx_map.items()]


        idxs.sort(key=lambda x: x[1])

        for acc, idx in idxs:  # Iterate through
            uparc_annots = uparc_xml_dict(uparc_annots_from_file(raw_uniparc, idx=idx))  # uparc annot dictionary
            db_refs = uparc_annots["dbReference"]  # List of databse entries mapped to this uparc

            # # TODO: For testing, DELETE
            # print("DB_REFS:")
            # print(db_refs)
            #
            # # TODO: For testing, DELETE
            # print("PRINTING REFS:")
            # for ref in db_refs:
            #     print(ref)

            if isinstance(db_refs, list):  # Multiple database references
                uprot_accs = [ref["@id"] for ref in db_refs if
                              (ref["@type"] == "UniProtKB/TrEMBL" or ref["@type"] == "UniProtKB/Swiss-Prot")
                               and ref["@active"] == 'Y']
            elif (db_refs["@type"] == "UniProtKB/TrEMBL" or db_refs["@type"] == "UniProtKB/Swiss-Prot") \
                    and db_refs["@active"] == 'Y':    # Single database reference
                uprot_accs = [db_refs["@id"]]
            else:
                uprot_accs = []

            mappings[acc] = uprot_accs

        if accs and report_missing:
            print(f"The following IDs were not found in the provided index file: {', '.join(missing_ids)}")

    else:
        with open(raw_uniparc) as in_file:
            xml_str = ""
            acc = None
            while True:
                line = in_file.readline()

                if not line.strip():  # EOF
                    break

                if line.startswith("<entry dataset"):   # Beginning of record
                    xml_str = line
                    acc_line = in_file.readline()   # Read accession line
                    acc = acc_line.strip().split('>')[1].split('<')[0]
                    xml_str += acc_line

                elif line.strip().startswith("</entry>"):  # End of record
                    if not accs or acc in accs:  # If acc in specified list OR we're taking all records
                        xml_str += line
                        uparc_annots = uparc_xml_dict(xml_str)   # Parse XML to uparc annot dictionary

                        db_refs = uparc_annots["dbReference"]
                        if isinstance(db_refs, list):  # Multiple database references
                            uprot_accs = [ref["@id"] for ref in db_refs if
                                          (ref["@type"] == "UniProtKB/TrEMBL" or ref["@type"] == "UniProtKB/Swiss-Prot")
                                          and ref["@active"] == 'Y']
                        elif (db_refs["@type"] == "UniProtKB/TrEMBL" or db_refs["@type"] == "UniProtKB/Swiss-Prot") and \
                                db_refs["@active"] == 'Y':  # Single database reference
                            uprot_accs = [db_refs["@id"]]
                        else:
                            uprot_accs = []

                        mappings[acc] = uprot_accs

                else:
                    xml_str += line

    # Write uparc2uprot file if specified
    if parc2prot_file:
        with open(parc2prot_file, 'w') as out_file:
            for uparc, uprots in mappings.items():
                out_file.write(f"{uparc}\t{','.join(uprots)}\n")

    # Write uprot2uparc file if specified
    if prot2parc_file:
        reverse_mappings = {}
        for uparc, uprots in mappings.items():
            for uprot in uprots:
                reverse_mappings[uprot] = uparc
        with open(prot2parc_file, 'w') as out_file:
            for uprot, uparc in reverse_mappings.items():
                out_file.write(f"{uprot}\t{uparc}\n")


def uprot_annots_from_uparc(
        uparc_annot_file,
        uprot_annot_file,
        uparc_idx=None,
        prot2parc=None,
        max_batch_size=200,
        max_threads=10
):
    """ Retrieve raw Uniprot annotations corresponding to all Uniparc annotations from a given raw annotation file.
    This function can also be used to generate index files and annotation maps for Uniparc annotations if they are not
    already available. If these files are already present and follow naming conventions, they need not be specified. """

    with open("test_debug.log", 'a') as log:
        log.write("Started uprot_annots_from_uparc\n")

    # Uniparc index file
    if not uparc_idx:  # May be present with default name or we may need to create it
        uparc_idx_name = uparc_annot_file.split('.txt')[0]+'.idx'
        if not os.path.exists(uparc_idx_name):  # Need to create
            uniparc_annot_idx(uparc_annot_file, uparc_idx_name)
    else:
        uparc_idx_name = uparc_idx

    # ID mapping
    if not prot2parc:
        parc2prot_name = uparc_annot_file.split('.txt')[0] + '_parc2prot.tsv'
        prot2parc_name = uparc_annot_file.split('.txt')[0] + '_prot2parc.tsv'
        if not os.path.exists(prot2parc_name):
            up_id_mapping(uparc_annot_file, parc2prot_file=parc2prot_name, prot2parc_file=prot2parc_name,
                          idx_file=uparc_idx_name)
    else:
        prot2parc_name = prot2parc

    # Fetch Uniprot annots
    fetch_raw_to_file(
        list(file_util.tab_del_file_to_dict(prot2parc_name).keys()),
        uprot_annot_file,
        failure_log=uprot_annot_file.split('.txt')[0]+'.log',
        max_batch_size=max_batch_size,
        max_threads=max_threads
    )

    # Create index file
    uniprot_annot_idx(
        uprot_annot_file,
        uprot_annot_file.split('.txt')[0]+'.idx'
    )


def retrieve_all_annots(
        uparc_seqs,
        uparc_annot_file,
        uprot_annot_file,
        max_threads=10
):
    """ Retrieve Uniparc annots and all annotations from available member databases (e.g. UniprotKB) which map to any
     Uniparc entries. Index and mapping files will be created following naming conventions of this module. Uniparc
     sequences can be provided as a list of identifiers or a fasta file from which IDs are extracted. """

    if not isinstance(uparc_seqs, list):  # File provided so must extract
        uparc_seqs = list([seq.name for seq in SeqIO.parse(uparc_seqs, "fasta")])

    # Fetch uparc annots
    fetch_raw_to_file(uparc_seqs, uparc_annot_file, dbName='uniparc', format='uniprotrdfxml',
                      failure_log=uparc_annot_file.split('.txt')[0]+'.log', max_threads=max_threads)

    # TODO: For debugging - DELETE
    with open("test_debug.log", 'a') as log:
        log.write("Finished fetch_raw_to_file\n")

    # Generate uparc index and map files + fetch mapped uprot sequences
    uprot_annots_from_uparc(uparc_annot_file, uprot_annot_file, max_threads=max_threads)


def fetch_annot(
        entryId,
        dbName='uniprotkb',
        format='uniprot',
        id_fix=True
):
    """ Fetch annotations for a single query ID for a given dbfetch database (see EBI website). Note, only specific
    database searches are handled here. See fetch_annots_raw for fetching non-processed information in text format for
    other dbfetch databases not yet handled by this function.
     Returns None if search fails (i.e. not a valid
     ID for the given database. Function adapted from binfpy.webservice """

    # Base url for all EBI web services
    __ebiUrl__ = 'http://www.ebi.ac.uk/Tools/'

    # Construct URL
    # Note: hard-coded to return raw text - can change to allow html in future
    url = __ebiUrl__ + 'dbfetch/dbfetch?style=raw&Retrieve=Retrieve&db=' + dbName + '&format=' + format + '&id=' + entryId


    # Get the entry
    # TODO: Currently just working with uniprotkb. Update for structural databses + Uniref + Uniparc
    try:
        data = urllib.request.urlopen(url).read().decode("utf-8")
        if data.startswith("ERROR"):
            return None  # Assume provided ID is not uniprotkb, so return None
        return data.split('\n')
    except urllib.error.HTTPError as ex:  # Assume non-query related issue (e.g. connection issue). Warn user but keep going
        raise RuntimeWarning(ex.read())


def fetch_annots(seqs, annots='all'):
    """ Fetch Uniprot annotations for a given set of Uniprot sequences. Input may be in the form of sequence names, a
     fasta file, or a list of Bio.SeqRecords. Returns a dictionary mapping sequence ID : Uniprot annotations.
     Annotations are provided for each of the annot types specified, corresponding to those defined in ANNOTS. If annots
     are None, returns the un-processed annotation string. If annots are 'all', all annotations defined in ANNOTS are
     used. If an error is encountered processing a given annotation, None is returned for that annotation.
     If search fails for a given sequence, it is mapped to None rather than a dictionary.
     TODO: Currently assumes all entries are Uniprotkb and returns None if search fails.
     Future additions to allow checking of ID type and subsequently querying the correct database."""

    # Get seq names
    if isinstance(seqs, list):
        if isinstance(seqs[0], SeqIO.SeqRecord):  # Need names from SeqRecords
            seqs = [seq.name for seq in seqs]
        # Else assume we already have seqnames
    else:  # Assume we have a Fasta file
        seqs = [seq.name for seq in list(SeqIO.parse(seqs, 'fasta'))]

    all_annots = {}

    for seq in seqs:
        raw_annots = fetch_annot(seq)    # TODO: Currently just trying UniprotKB and returning none otherwise

        if raw_annots:  # We got annot data

            if not annots:  # Just want the raw annots
                all_annots[seq] = raw_annots

            else:

                this_annots = {}

                if annots == 'all':   # Want all annots defined in UP_ANNOTS
                    annots = list(UP_ANNOTS.keys())

                for annot in annots:  # Apply all specified annots
                    try:
                        this_annots[annot] = UP_ANNOTS[annot](raw_annots)
                    except:
                        this_annots[annot] = None

                all_annots[seq] = this_annots

        else:   # Couldn't get annots so just put None
            all_annots[seq] = None

    return all_annots


def fetch_annots_to_file(out_file, seqs, annots='all'):
    """ Fetch specified annots and write to Figree-formatted annotation file. Missing annotations are written as a space
    """
    annot_dict = fetch_annots(seqs, annots)
    annot_file_from_dict(out_file, annot_dict)


def annot_file_from_dict(out_file, annot_dict):
    """ Output an annotation file in Figtree (tab-delimited) format from a nested annotations dictionary as output by
     fetch_annots. I.e. annot_dict format is {seq_1 : {annot_1 : value_1, ... }, ... }. Where a given annotation is not
     available for a sequence, include a space character in the annotation file (otherwise Figree doesn't like it)."""

    with open(out_file, 'w') as file:

        # Collect all annotation types
        annot_types = set()
        for annots in annot_dict.values():
            if not annots:  # No annots for this seq
                continue
            annot_types.update(annots.keys())
        annot_types = list(annot_types)
        annot_types.sort()

        # Write header
        file.write('taxa')
        for annot in annot_types:
            file.write(f'\t{annot}')
        file.write('\n')

        for seq, annots in annot_dict.items():

            if not annots:
                annots = {annot : ' ' for annot in annot_types}

            file.write(seq)

            for annot in annot_types:
                try:
                    annot_val = annots[annot]
                except KeyError:
                    annot_val = ' '
                except AttributeError:
                    annot_val = ' '
                if not annot_val:
                    annot_val = ' '
                file.write(f'\t{annot_val}')

            file.write('\n')


def annot_dict_from_file(annot_file):
    """ Return annotations from a given file in the standard nested dictionary format. """

    annot_dict = {}

    annots = [line.strip('\n') for line in open(annot_file).readlines()]  # Just split on \n as there may be spaces
    annot_types = annots[0].split('\t')[1:]  # First col is just header for seq names

    for seq_annots in annots[1:]:

        vals = seq_annots.split('\t')
        name = vals[0]

        # Make empty values (single space) from file = None
        for i in range(len(vals)):
            if vals[i] == ' ':
                vals[i] = None

        annot_dict[name] = {}

        if len(vals[1:]) != len(annot_types):
            raise RuntimeError(f'Incorrect number of annotations for sequence: {name}.')

        for i in range(len(annot_types)):
            annot_dict[name][annot_types[i]] = vals[i+1]  # Add one to index becas

    return annot_dict


def clean_annot_dict(annot_dict):
    """ Clean annotation dictionary by replacing empty or missing values as None."""

    # Get union of all annotation types
    annot_types = set()
    for annots in annot_dict.values():
        if not annots:  # This seq doesn't have annots
            continue
        annot_types.update(annots.keys())
    annot_types = list(annot_types)
    annot_types.sort()

    all_cleaned = {}
    for seq, annots in annot_dict.items():
        this_cleaned = {}
        for annot in annot_types:
            try:
                annot_val = annots[annot]
            except KeyError:
                annot_val = None
            except TypeError:
                annot_val = None
            except AttributeError:
                annot_val = None
            if annot_val == ' ':  # Handles blank space only annotations in file
                annot_val = None
            this_cleaned[annot] = annot_val
        all_cleaned[seq] = this_cleaned

    return all_cleaned


def merge_annot_dicts(annot_dicts):
    """ Merge annotation dictionaries of the form output by above functions. Blank annots are converted to None. Missing
     annotation values for a given sequence are imputed as None. For contradicting annotation values, priority is
     given to the first listed annot_dict. """

    if isinstance(annot_dicts, dict):  # Just one dict
        return annot_dicts

    merged = annot_dicts[-1]  # Start at lowest priority, iteratively merging higher priority sets in

    for annot_set in annot_dicts[-2::-1]:
        annot_set = clean_annot_dict(annot_set)  # Make sure every annot has None values for missing annots
        for seq, annots in annot_set.items():
            if seq not in merged.keys():  # New sequence
                merged[seq] = annots
            else:
                for annot_type, annot_val in annots.items():
                    if annot_val:   # Don't overwrite existing annotations with None
                        merged[seq][annot_type] = annot_val

    merged = clean_annot_dict(merged)

    return merged


def merge_annots(
        annot_dicts=None,
        annot_files=None,
        out_file=None,
        no_return=True
):
    """ Combine annotations from multiple file and/or dictionary sources. Merged annotations may be written out to
    a new annotations file, returned as a merged dictionary, or both.
    NOTE: Priority for merging contradicting annotations is: dicts first to last index, followed by files first to last
    index. If an existing file should take precedence, read it in as a dictionary first, and assign it appropriate
    priority. """

    if not annot_dicts:
        annot_dicts = []
    if not annot_files:
        annot_files = []

    if not annot_files and not annot_dicts:
        raise RuntimeError("No annotations provided.")

    if no_return and not out_file:
        raise RuntimeError("If no return is selected, an output file name must be specified.")

    if isinstance(annot_dicts, dict):
        annot_dicts = [annot_dicts]

    if isinstance(annot_files, str):
        annot_files = [annot_files]

    # Read in annot files as dicts
    for file in annot_files:
        annot_dicts.append(annot_dict_from_file(file))

    merged = merge_annot_dicts(annot_dicts)

    if out_file:
        annot_file_from_dict(out_file, merged)

    if not no_return:
        return merged


def create_itol_metadata(
        itol_file,
        annot_file=None,
        annot_dict=None,
        exclude_fields=None
):
    """ Write an ITOL metadata file from a tab-delimited annotation file or annotation dictionary. Selected annotation
     fields can be excluded. """

    if not (annot_file or annot_dict):
        raise RuntimeError("Either an annotation file or dictionary must be provided.")

    if not exclude_fields:
        exclude_fields = ["sequence"]  # Always want to exlclude sequence

    if annot_file:  # Read in file if provided
        annot_dict = annot_dict_from_file(annot_file)

    annot_dict = clean_annot_dict(annot_dict)  # Make all missing/empty values None

    # Get all annotation fields
    for node, annots in annot_dict.items():
        fields = list(annots.keys())
        break

    # Exclude specific fields as desired
    for field in exclude_fields:
        if field in fields:
            fields.remove(field)

    # See ITOL website for metadata file format
    with open(itol_file, 'w') as out_file:
        out_file.write("METADATA\n")
        out_file.write("SEPARATOR TAB\n")
        out_file.write("FIELD_LABELS\t" + '\t'.join(fields) + '\n')
        out_file.write("DATA\n")
        for node, annots in annot_dict.items():
            # Null values written as .
            out_file.write(f"{node}\t" + '\t'.join([annots[field] if annots[field] else '.' for field in fields]) + '\n')


def itol_internal_labels(
        file_name,
        nodes,
        color="#0000ff",
        legend_name="Ancestors"
):
    """ Write an itol dataset_symbol formatted file which places shapes on internal nodes of interest. """

    if not isinstance(nodes, list):
        nodes = [nodes]

    with open(file_name, 'w') as file:

        file.write(f"DATASET_SYMBOL\nSEPARATOR TAB\nDATASET_LABEL\t{legend_name}\nCOLOR\t{color}\nMAXIMUM_SIZE\t10\nDATA")
        for node in nodes:
            file.write(f"\n{node}\t2\t8\t{color}\t1\t1")


def itol_strip_from_annot(
        full_annot_file,
        color_strip_file,
        annot_label,
        binary_annot=False,
        binary_color="#ff0000"
):
    """ Create an iTol colour strip annotation file using a data for a particular annotation label in an existing
    annotation (.annot) file. If a binary annotation is specified, specific annotation values are ignored and any
    sequence with a non-null annotation are tagged with the same colour and with label 'True'. """

    # Colour scheme for annotation label groups
    tab20 = mpl.colormaps["tab20"]
    hex_cols = [mpl.colors.to_hex(c) for c in tab20.colors]

    # Sort annotations for given label
    full_annots = annot_dict_from_file(full_annot_file)
    spec_annots = {seq : seq_annots[annot_label]
                   for seq, seq_annots in full_annots.items()}

    if binary_annot:

        spec_annots = {seq : "True"
                       for seq, val in spec_annots.items() if val}

        color_map = {"True" : binary_color}

    else:  # Allocate specific colours to unique annot values

        # Rank specific values by frequency
        val_counts = {}
        for val in spec_annots.values():
            try:
                val_counts[val] += 1
            except KeyError:
                val_counts[val] = 1
        val_counts = list(val_counts.items())
        val_counts.sort(key=lambda x:x[1], reverse=True)

        # Give specific labels/colours to top 10, then collapse into "Other"
        color_map = {}
        top_10_vals = []

        for i in range(min(10, len(val_counts))):
            color_map[val_counts[i][0]] = hex_cols[i]
            top_10_vals.append(val_counts[i][0])

        if len(val_counts) > 10:
            # Need to replace all other values with "Other"
            for seq in spec_annots:
                if spec_annots[seq] not in top_10_vals:
                    spec_annots[seq] = "Other"
            color_map["Other"] = hex_cols[10]

    # Write annotation file
    with open(color_strip_file, 'w') as out_f:

        out_f.write(f"DATASET_COLORSTRIP\nSEPARATOR COMMA\n")
        out_f.write(f"DATASET_LABEL,{annot_label}\nCOLOR,#ff0000\n")

        # The actual data
        out_f.write("DATA\n")
        for seq, annot_val in spec_annots.items():
            out_f.write(f"{seq},{color_map[annot_val]},{annot_val}\n")


def custom_db_annots(seqs, function):
    """ Return a dict mapping seq_name to annotation value as determined by a lambda function. Input may be in the form
    of sequence names, a fasta file, or a list of Bio.SeqRecords. Lambda function should take a single input
    to which the Uniprot annotation string for each sequence (if available) will be supplied. """

    if isinstance(seqs, list):
        if isinstance(seqs[0], SeqIO.SeqRecord):   # Need to convert to names
            seqs = [seq.name for seq in seqs]
        # Else assume that they're already accessions

    else:  # Assume seqs is fata file name
        seqs = [seq.name for seq in list(SeqIO.parse(seqs, 'fasta'))]

    full_annots = fetch_annots(seqs)

    target_annots = {}
    for name, annots in full_annots.items():
        try:
            target_annots[name] = function(annots)
        except:
            target_annots[name] = ' '  # Blank where given annotation isn't available

    return target_annots


def uprot_annots_from_file(
        raw_annots,
        idx=None,
        seq_id=None,
        id_type='AC'
):
    """ Extract a single Uniprot annotation record from an existing file, using offset index if available. Returns
     None if record not found for given accession. Note: very sensitive to any changes to raw annotation file
     structure. """

    if not (idx or seq_id):
        raise RuntimeError("Either an index or sequence ID must be provided for annotation retrieval. ")

    with open(raw_annots) as annots:

        if idx:
            annots.seek(idx)
            txt = ""
            while True:  # Wait for a line with just '//' representing end of record
                line = ""
                while True:  # Append record line by line
                    char = annots.read(1)
                    line += char
                    if char == '\n':
                        txt += line
                        break
                if line == "//\n":
                    break
            return txt

        else:
            txt = ""
            found_seq = False
            while True:
                line = annots.readline()
                if not line.strip():  # At EOF, seq not found
                    return None
                elif line[:2] == id_type:  # Prefix of annotation line to match to
                    if seq_id in line:  # This is the entry we want
                        found_seq = True
                if line.strip() == "//":  # End of record
                    if found_seq:
                        return txt
                    else:
                        txt = ""
                        continue
                txt += line


def uparc_annots_from_file(raw_annots, idx=None, seq_id=None):
    """ Extract the raw XML text for a single Uniparc sequence as a string from an annotation file. """

    if not (idx or seq_id):
        raise RuntimeError("Either an index or Uniparc sequence accession must be provided.")

    with open(raw_annots) as annots:

        if idx:
            annots.seek(idx)
            txt = ""
            while True:  # Wait for '</rdf:RDF>' indicating end of entry
                line = annots.readline()
                txt += line
                if line.strip().startswith("</rdf:RDF>"):
                    break
            return txt

        else:
            txt = ""
            while True:
                line = annots.readline()
                if not line.strip():  # At EOF, seq not found
                    return None
                elif line.strip().startswith("<rdf:RDF"):   # Start of record
                    txt = ""
                elif line.strip().startswith("</rdf:RDF>"):  # End of record
                    txt += line
                    if f'rdf:about="{seq_id}"' in txt:
                        return txt
                txt += line


def uparc_xml_dict(xml_str):
    """ From raw XML string, extract a Uniparc annotation dictionary. Returns annotations as an Ordered Dictionary. """

    return xmltodict.parse(xml_str)["entry"]


def extract_uparc_sources(xml_str, db_map_file=None):
    """ From a uniprotrdfxml string, extract all accessions from of source databases (note: not all databases currently
    supported) and return as a dictionary as {database_name : [acc1, acc2], .....}. """

    # Note: listed dbs are those which follow the standard format. Any with exceptions are handled specifically below
    supported_dbs = ["uniprot"]

    # TODO: For testing, remove try/except
    try:
        xml_dict = xmltodict.parse(xml_str)
    except:
        print(xml_str)
        raise RuntimeError(f"Failed")

    sources = xml_dict['rdf:RDF']["rdf:Description"][0]["sequenceFor"]  # Get list of source URLs
    if not isinstance(sources, list):
        sources = [sources]

    source_dict = {}

    for source in sources:
        try:
            url = source["@rdf:resource"]

            if "/pdb/" in url:  # Different url format to other DBs
                if "http" in url and url.split("/")[-4] == "pdb":
                    if "pdb" in source_dict.keys():
                        source_dict["pdb"].append(url.split("/")[-3])
                    else:
                        source_dict["pdb"] = [url.split("/")[-3]]

            else:
                if "http" in url and url.split("/")[-2] in supported_dbs:
                    db = url.split("/")[-2]
                    if db in source_dict.keys():
                        source_dict[db].append(url.split("/")[-1])
                    else:
                        source_dict[db] = [url.split("/")[-1]]

        except IndexError:
            print(f"Error extracting source entry: {url}")

        return source_dict if source_dict else None


def extract_uparc_proteomes(raw_uparc_annots, select_ids=None):
    """ From raw Uniparc annotations, extract all mapped Uniprot proteome accessions (UPO...) and return as a set. """

    proteome_dict = {}

    # Create index file if not already available
    idx_file = raw_uparc_annots.split(".")[0] + ".idx"
    if idx_file not in os.listdir():
        uniparc_annot_idx(raw_uparc_annots, idx_file)

    idxs = file_util.tab_del_file_to_dict(idx_file)

    if not select_ids:  # If specific seqs aren't specified, proceed for all in the annotation file
        select_ids = list(idxs.keys())

    for seq_id, idx in idxs.items():

        if seq_id in select_ids:

            data = xmltodict.parse(uparc_annots_from_file(raw_uparc_annots, idx=int(idx)))["rdf:RDF"]["rdf:Description"]

            # Check which data sources are associated with a Uniprot proteome
            for data_dict in data:
                if "proteome" in data_dict:
                    proteome_id = data_dict["proteome"]["@rdf:resource"].split("/proteomes/")[1].split("#")[0]
                    try:
                        proteome_dict[seq_id].add(proteome_id)
                    except KeyError:
                        proteome_dict[seq_id] = {proteome_id}

    return proteome_dict


def assembly_from_proteome(up_proteome_ids, batch_size=100):
    """ From a list of Uniparc proteomes, fetch associated genome assemblies. """

    url = "https://rest.uniprot.org/proteomes/search"

    results = {}

    for i in range(0, len(up_proteome_ids), batch_size):
        chunk = up_proteome_ids[i:i + batch_size]
        query = " OR ".join(f"upid:{u}" for u in chunk)

        params = {
            "query": query,
            "format": "json",
            "fields": "upid,genome_assembly",
            "size": 500
        }

        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()

        for entry in data.get("results", []):
            results[entry["id"]] = entry.get("genomeAssembly", {}).get("assemblyId")

    return results


def gtdb_metadata_from_assemblies(assembly_ids, gtdb_metadata_file=None):
    """ Extract GTDB metadata mapped to an NCBI genome assembly ID. If the genome is not present, it is skipped with
     no error raised. """

    if not gtdb_metadata_file:  # Assume file has default name
        gtdb_metadata_file = "bac120_metadata.tsv"

    # Source GTDB metadata file if it isn't present
    if gtdb_metadata_file not in os.listdir():

        if gtdb_metadata_file + ".gz" in os.listdir():  # Leave decompressed for future use
            with gzip.open(gtdb_metadata_file + ".gz", "rb") as in_file:
                with open(gtdb_metadata_file, "wb") as out_file:
                    shutil.copyfileobj(in_file, out_file)

        else:  # Download metadata
            gtdb_metadata_file = "bac120_metadata.tsv"
            url = "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata_tsv.gz"

            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open("bac120_metadata.tsv.gz", "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

            # Decompress
            with gzip.open(gtdb_metadata_file + ".gz", "rb") as in_file:
                with open(gtdb_metadata_file, "wb") as out_file:
                    shutil.copyfileobj(in_file, out_file)

    fields = [ # TODO: Add additional relevant fields for extraction and/or allow custom selection of fields
        "ncbi_genbank_assembly_accession",
        "ncbi_isolation_source"
    ]

    # Note: GTDB does not resolve at kingdom level
    tax_levels = {"d" : "domain", "p" : "phylum", "c" : "class", "o" : "order", "f" : "family",
                  "g" : "genus", "s" : "species"}

    results = {}

    with open(gtdb_metadata_file) as in_file:
        header = True
        for line in in_file:

            data = line.strip().split("\t")

            if header:
                n_cols = len(data)
                gtdb_id_idx = data.index("accession")
                assembly_id_idx = data.index("ncbi_genbank_assembly_accession")
                gtdb_tax_idx = data.index("gtdb_taxonomy")
                field_idxs = [data.index(field) for field in fields]
                header = False

            elif len(data) == n_cols and data[assembly_id_idx] in assembly_ids:  # Only extract annots for selected assemblies

                assembly = data[assembly_id_idx]

                # GTDB accession
                gtdb_acc = data[gtdb_id_idx]
                try:
                    results[assembly]["gtdb_accession"] = gtdb_acc
                except:
                    results[assembly] = {"gtdb_accession" : gtdb_acc}

                # Split GTDB taxonomy
                tax_str = data[gtdb_tax_idx]
                tax_cats = tax_str.split(";")
                for cat in tax_cats:
                    level, label = tax_levels[cat.split("__")[0]], cat.split("__")[1]
                    try:
                        results[assembly][f"gtdb_{level}"] = label
                    except KeyError:
                        results[assembly] = {f"gtdb_{level}" : label}

                # Extract other fields
                for i in range(len(field_idxs)):

                    try:
                        results[assembly][fields[i]] = data[field_idxs[i]]
                    except KeyError:
                        results[assembly] = {fields[i] : data[field_idxs[i]]}

    return results


def uparc_to_gtdb_annots(
        raw_uparc_annots,
        select_ids=None,
        gtdb_metadata_file=None
):
    """ For specified Uniparc annotations, map to GTDB metadata via assembly IDs (where available) and return a
    annotation dictionary in the standard internal format of this module. Note: the annotations fetched are currently
    dictated by the function directly extracting from the GTDB metadata tsv. """

    uparc_idx = raw_uparc_annots.split(".")[0] + ".idx"
    if uparc_idx not in os.listdir():
        uniparc_annot_idx(raw_uparc_annots, uparc_idx)

    # Get sets of Uniprot proteomes mapping to each Uniparc ID
    proteome_mappings = extract_uparc_proteomes(raw_uparc_annots, select_ids=select_ids)

    # Get genome assemblies associated with each mapped proteome
    # TODO: Currently taking the first proteome only for Uniparc IDs mapped to multiple proteomes
    proteome_mappings_single = {uparc_id : list(proteomes)[0] for uparc_id, proteomes in proteome_mappings.items()}
    assembly_mappings = assembly_from_proteome(list(proteome_mappings_single.values()))

    # Bridge Uniparc to assembly mappings
    uparc_to_assembly = {}
    for up_id, proteome_id in proteome_mappings_single.items():
        if proteome_id in assembly_mappings:
            uparc_to_assembly[up_id] = assembly_mappings[proteome_id]

    # Extract GTDB annotations for all available assemblies and map to relevant Uniparc ID
    gtdb_annots = gtdb_metadata_from_assemblies(list(uparc_to_assembly.values()))
    uparc_to_gtdb = {}
    for up_id, assembly_id in uparc_to_assembly.items():
        if assembly_id in gtdb_annots:
            uparc_to_gtdb[up_id] = gtdb_annots[assembly_id]

    return uparc_to_gtdb


def tax_proportion_from_annots(
        seq_ids,
        seq_annots,
        hc_map=None,
        hc_name=None,
        hc_level=None,
        weight_by_cluster_size=True,
        ncbi_tax_prefix="tax",
        gtdb_tax_prefix="gtdb"
):
    """ Perform taxonomic proportion analysis for a given set of sequences. Attains both NCBI and GTDB taxonomy
     proportions. If a hierarchical clustering instance is specified, proportions are calculated on the basis of all
    sequences clustered with the representative, with optional down-weighting for cluster size if desired to account
    for biases in taxon sampling. """

    if hc_name and not hc_map:

        if not hc_level:
            raise RuntimeError("The hierarchical clustering level must be specified.")

        hc_map = cluster.get_hc_maps(name=hc_name)[hc_level]

        # Check that all specified sequences are reps at the chosen hc level
        for seq in seq_ids:
            if seq not in hc_map:
                raise RuntimeError(f"Sequence {seq} is not a representative at the specified clustering threshold.")

    if isinstance(seq_annots, str):  # Assume annots are in a .annot file
        seq_annots = annot_dict_from_file(seq_annots)
    elif not isinstance(seq_annots, dict):
        raise RuntimeError("Annotations must be provided either as a .annot file or a properly formatted nested "
                           "annotation dictionary.")

    tax_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    tax_prefix = [ncbi_tax_prefix, gtdb_tax_prefix]
    results = [{level : {} for level in tax_levels} for i in range(len(tax_prefix))]  # NCBI, GTDB
    total_cluster_members = 0
    non_gtdb_seqs = 0
    non_gtdb_clusters = 0

    for rep in seq_ids:

        # Per-cluster label counts for NCBI and GTDB in each taxonomy level
        cluster_results = [{level : {} for level in tax_levels} for i in range(len(tax_prefix))]
        found_gtdb = False

        if hc_map:  # Need to analyse per cluster
            cluster_seqs = hc_map[rep]
        else:  # Treat each rep as a single member cluster
            cluster_seqs = [rep]

        total_cluster_members += len(cluster_seqs)

        for seq in cluster_seqs:

            if seq in seq_annots:

                # TODO: Change to checking presence of gtdb_accession for cluster mapping instead of ncbi
                if seq_annots[seq]["ncbi_genbank_assembly_accession"]:
                    found_gtdb = True
                else:
                    non_gtdb_seqs += 1

                for tax_level in tax_levels:

                    # NCBI taxonomy, then GTDB taxonomy
                    for i in range(len(tax_prefix)):
                        label = seq_annots[seq][f"{tax_prefix[i]}_{tax_level}"]
                        if label:
                            try:
                                cluster_results[i][tax_level][label] += 1
                            except KeyError:
                                cluster_results[i][tax_level][label] = 1

        if not found_gtdb:
            non_gtdb_clusters += 1

        # Add cluster tax counts to overall counts
        for tax_level in tax_levels:
            for i in range(len(tax_prefix)):
                if cluster_results[i][tax_level]:  # There may not be any annotations at this level - if so, ignore

                    # Apply normalisation if weight_by_cluster_size specified
                    norm_factor = sum([cnt for cnt in cluster_results[i][tax_level].values()]) \
                        if weight_by_cluster_size else 1

                    for label, cnt in cluster_results[i][tax_level].items():

                        try:
                            results[i][tax_level][label] += cnt / norm_factor
                        except KeyError:
                            results[i][tax_level][label] = cnt / norm_factor

        # Calculate proportions
        tax_props = [{level : {} for level in tax_levels} for i in range(len(tax_prefix))]
        for tax_level in tax_levels:
            for i in range(len(tax_prefix)):
                if results[i][tax_level]:
                    total_cnts = sum([cnt for cnt in results[i][tax_level].values()])
                    for label, cnt in results[i][tax_level].items():
                        tax_props[i][tax_level][label] = cnt / total_cnts

    return tax_props, (len(seq_ids), total_cluster_members), (non_gtdb_clusters, non_gtdb_seqs)


def clade_tax_comparison(
        t,
        seq_annots,
        clade_seqs,
        tax_rank,
        out_seqs=None,
        clade_names=None,
        hc_map=None,
        hc_name=None,
        hc_level=None,
        weight_by_cluster_size=True,
        ncbi_tax_prefix="tax",
        gtdb_tax_prefix="gtdb",
        itol_piechart_prefix=None,
        itol_dataset_label=None
):
    """ Perform taxonomic distribution comparison between NCBI and GTDB taxonomy for specified clades in a phylogeny.
     Clades should be specified as o"""

    rank_props = []  # Indexed by order of clade_seqs and out_seqs

    t = tree.load_tree(t)

    if not out_seqs:
        out_seqs = [None for _ in range(len(clade_seqs))]

    if not clade_names:
        clade_names = [None for _ in range(len(clade_seqs))]

    for i in range(len(clade_seqs)):

        if out_seqs[i]:
            t = tree.outgroup_root(t, clade_seqs[i], out_seqs[i])

        sub_seqs = tree.get_subtree_leaves(t, ext_nodes=clade_seqs[i])

        tax_prop_results = tax_proportion_from_annots(sub_seqs, seq_annots, hc_map=hc_map, hc_name=hc_name,
                                                  hc_level=hc_level, weight_by_cluster_size=weight_by_cluster_size,
                                                  ncbi_tax_prefix=ncbi_tax_prefix, gtdb_tax_prefix=gtdb_tax_prefix)

        print(f"{clade_names[i] if clade_names[i] else 'Clade '+str(i+1)}")
        print(f"{tax_prop_results[2][1]} of {tax_prop_results[1][1]} total sequences not mapped to a GTDB genome.")
        print(f"{tax_prop_results[2][0]} of {tax_prop_results[1][0]} total clusters not mapped to any GTDB genomes.")
        print()

        rank_props.append(
            ([tax_prop_results[0][0][tax_rank], tax_prop_results[0][1][tax_rank]],   # [NCBI, GTDB]
            tax_prop_results[1],
            tax_prop_results[2]
                           ))

    # Write itol piechart annotation files if specified
    if itol_piechart_prefix:

        tax_types = [ncbi_tax_prefix, gtdb_tax_prefix]

        for tax_type_idx in range(len(tax_types)):

            if not itol_dataset_label:
                itol_dataset_label = f"{tax_types[tax_type_idx]}_{tax_rank}"

            # Colour scheme for charts and legend
            tab20 = mpl.colormaps["tab20"]
            hex_cols = [mpl.colors.to_hex(c) for c in tab20.colors]

            # Order tax labels in legend by overall abundance across all clades
            total_abundance = {}
            for j in range(len(clade_seqs)):
                # print(f"{j}\t{len(clade_seqs)}\t{len(rank_props)}")
                # print(rank_props[j])
                # print(rank_props[j][0])
                # print(rank_props[j][0][tax_type_idx])
                for label, prop in rank_props[j][0][tax_type_idx].items():
                    try:
                        total_abundance[label] += prop
                    except KeyError:
                        total_abundance[label] = prop
            total_abundance = list(total_abundance.items())
            total_abundance.sort(reverse=True, key=lambda x: x[1])
            # If there are more than 10 labels total, collapse all but top10 into "Other"
            top10 = [total_abundance[i][0] for i in range(10 if len(total_abundance) > 10 else len(total_abundance))]

            # Proportions for pie chart, including anything outside top 10 collapsed into "Other"
            chart_props = []  # indexed by order of clade_seqs
            for j in range(len(clade_seqs)):
                chart_props.append({})
                for label, prop in rank_props[j][0][tax_type_idx].items():
                    if label in top10:
                        chart_props[j][label] = prop
                    else:
                        try:
                            chart_props[j]["Other"] += prop
                        except KeyError:
                            chart_props[j]["Other"] = prop

            # Write itol metadata file
            with open(f"{itol_piechart_prefix}_{tax_types[tax_type_idx]}.itol", 'w') as itol_file:

                itol_file.write(f"DATASET_PIECHART\nSEPARATOR COMMA\nDATASET_LABEL,{itol_dataset_label}\nCOLOR,#ff0000\n")

                if len(total_abundance) > 10:
                    n_cats = 11
                    chart_labels = top10 + ["Other"]
                else:
                    n_cats = len(total_abundance)
                    chart_labels = top10

                chart_colors = [hex_cols[j] for j in range(n_cats)]

                itol_file.write(f"FIELD_COLORS,{','.join(chart_colors)}\n")

                itol_file.write(f"FIELD_LABELS,{','.join(chart_labels)}\n")

                # Legend
                itol_file.write(f"LEGEND_TITLE,{tax_rank}\nLEGEND_SCALE,1\nLEGEND_SHAPES,"
                                f"{','.join(['2' for _ in range(n_cats)])}\n")
                itol_file.write(f"LEGEND_COLORS,{','.join(chart_colors)}\n")
                itol_file.write(f"LEGEND_LABELS,{','.join(chart_labels)}\n")

                # Data
                itol_file.write("DATA\n")

                for j in range(len(clade_seqs)):
                    itol_file.write(f"{clade_seqs[j][0]}|{clade_seqs[j][1]},-1,80")
                    for label in top10:
                        if label in chart_props[j]:
                            itol_file.write(f",{str(chart_props[j][label])}")
                        else:
                            itol_file.write(f",0")
                    if "Other" in chart_props[j]:
                        itol_file.write(f",{str(chart_props[j]['Other'])}")
                    itol_file.write("\n")

    return rank_props


def extract_annots(
        raw_annots,
        seqs,
        functions,
        idx_file=None,
        record_type="uniprot",
        id_type='AC',
        out_file=None,
        no_return=True,
        field_prefix=None,
        max_threads=1,
        report_failed_retrieval=False
):
    """ Extract annotation records for several sequences from a raw annotation file using offset index if available.
    Run defined functions over each to attain custom annotations. If output filename is provided write annotations in
    tab-delimited format, else return as a dictionary. Prints warnings for sequences which are not found. Functions for
    custom annotations should be supplied as a dictionary mapping annotation label to a lambda function acting on a list
    of lines for raw Uniprot annotation records or an XML ordered dictionary for Uniparc annotations. """

    if no_return and not out_file:
        raise RuntimeError("Annotations must either be written to file or returned as a dictionary.")

    annot_dict = {}

    if idx_file:  # Read in index if available
        with open(idx_file) as file:
            idxs = {line.strip().split('\t')[0] : int(line.strip().split('\t')[1])
                    for line in file.readlines()}

    # Function for extracting annots for a single sequence
    def get_annots_for_seq(seq, lock):

        if record_type == 'uniprot':
            annot_txt = uprot_annots_from_file(
                raw_annots,
                idx=idxs[seq] if idx_file else None,
                seq_id=seq,
                id_type=id_type
            )

        elif record_type == 'uniparc':
            annot_txt = uparc_annots_from_file(
                raw_annots,
                idx=idxs[seq] if idx_file else None,
                seq_id=seq
            )

        else:
            raise RuntimeError(f"Annotation fetching for record type {record_type} not implemented.")

        if not annot_txt:  # If None retuned, sequence not found
            print(f"Sequence {seq} not found in annotation file. Specified ID type may be incorrect.")

        else:  # Annotation record was found
            annot_dict[seq] = {}

            # Process raw annotation text into format appropriate for downstream annotation extraction
            if record_type == 'uniprot':
                annot_split = annot_txt.splitlines()  # uprot lambda function input is list of lines
            elif record_type == 'uniparc':
                annot_split = uparc_xml_dict(annot_txt)  # uparc lambda function input is xml ordered dict

            for label, function in functions.items():
                try:
                    if field_prefix:  # Apply prefix to field label if one is provided
                        label = field_prefix + label
                    with lock:
                        annot_dict[seq][label] = function(annot_split)

                except:
                    if report_failed_retrieval:  # Print warning if annotation function fails for any reason
                        print(f"Annotation {label} retrieval failed for sequence {seq}.")
                    else:
                        pass

    # Set up for multi-threading
    executor = futures.ThreadPoolExecutor(max_workers=max_threads)
    lock = threading.Lock()
    future_store = [executor.submit(
        get_annots_for_seq,
        seq,
        lock
    ) for seq in seqs]
    futures.wait(future_store)

    annot_dict = clean_annot_dict(annot_dict)   # Make missing annots None

    if out_file:
        annot_file_from_dict(out_file, annot_dict)

    if not no_return:
        return annot_dict


def map_and_extract_annots(
        parc2prot_map,
        uprot_annot_file,
        functions,
        out_file=None,
        no_return=False,
        uprot_idx=None,
        uparc_ids=None,
        full_sprot_map_file=None,
        record_type='uniprot',
        field_prefix=None,
        sprot_prefix=None,
        max_threads=2):
    """ Retrieve mapped annotations corresponding to a selection of Uniparc IDs (or all in available file if not
    specified). Returns an annotation dictionary and/or writes an annotation file with Uniparc IDs mapped to annotations
     of associated entries from other databases. An index file for the member annotation file will be created if not
     available. If the file is available and obeys naming convention, it need not be specified.
      """

    if no_return and not out_file:
        raise RuntimeError("An output file must be specified if the annotation dictionary is not returned.")

    # Get map of desired Uniparc IDs to a single active member database ID (if one exists)
    parc2prot_dict = file_util.tab_del_file_to_dict(parc2prot_map)
    single_map = {}  # Map for only desired Uniparc IDs with associated Uniprot entries
    if not uparc_ids:
        uparc_ids = list(parc2prot_dict.keys())

    # Keep track of UParc IDs mapped to extracted SProt sequences
    extracted_sprot = []   # NOTE: actually Uniparc accessions

    for seq in uparc_ids:

        try:
            mapped_uprot_str = parc2prot_dict[seq]
        except KeyError:
            print(f"Uniparc ID {seq} not present in map file.")
            continue

        if mapped_uprot_str:  # Only if it maps to at least one uniprot entry

            # Get all Uprot IDs (SProt + trembl)
            mapped_uprot_ids = mapped_uprot_str.split(',')

            # If Uniparc <--> SProt mapping is available, prioritise SProt over trembl
            if full_sprot_map_file:
                if ".sqlite" in full_sprot_map_file:  # Need to read from .tsv , not .sqlite
                    full_sprot_map_file = full_sprot_map_file.rsplit('.', maxsplit=1)[0] + ".tsv"
                # TODO: Inefficient doing this for every mapped Uparc ID - collect Sprot IDs once in outer scope
                mapped_sprot_ids = []
                with open(full_sprot_map_file) as sprot_f:
                    for line in sprot_f:
                        this_sprot_id = line.strip().split('\t')[0]
                        if this_sprot_id in mapped_uprot_ids:
                            mapped_sprot_ids.append(this_sprot_id)
                if mapped_sprot_ids:
                    single_map[seq] = mapped_sprot_ids[0]
                    extracted_sprot.append(seq)
                else:  # No SProt IDs so just take first trembl ID
                    single_map[seq] = mapped_uprot_ids[0]

            else: # Only one mapped UniProt or no SProt mapping file
                single_map[seq] = mapped_uprot_ids[0]

    # Create Uniprot annotation index if not available already
    if not uprot_idx:
        uprot_idx_name = uprot_annot_file.split('.txt')[0]+'.idx'
        if not os.path.exists(uprot_idx_name):
            uniprot_annot_idx(uprot_annot_file, uprot_idx_name)
    else:
        uprot_idx_name = uprot_idx

    # Fetch Uniprot annotations as dictionary
    uprot_annot_dict = extract_annots(
        uprot_annot_file,
        list(single_map.values()),
        functions,
        idx_file=uprot_idx_name,
        record_type=record_type,
        no_return=False,
        field_prefix=field_prefix,
        max_threads=max_threads
    )

    # Swap in Uniparc IDs for the member entry they map to
    uparc_annot_dict = {uparc_id : uprot_annot_dict[uprot_id]
                        for uparc_id, uprot_id in single_map.items()}

    # Tag Swissprot annotations explicitly
    if sprot_prefix:
        sprot_prefix = sprot_prefix.rsplit("_", maxsplit=1)[0] + "_"
        for seq in extracted_sprot:
            tagged_annots = {}
            for label, annot_val in uparc_annot_dict[seq].items():
                tagged_annots[sprot_prefix + label] = annot_val
            uparc_annot_dict[seq].update(tagged_annots)

    if out_file:
        annot_file_from_dict(out_file, uparc_annot_dict)
    elif not no_return:
        return uparc_annot_dict


def extract_sprot_ids(uprot_annot_file):
    """ Extract all Uniprot accessions corresponding to a Swissprot entry as a list from a raw annotation file. """

    sprot_ids = set()

    with open(uprot_annot_file) as file:
        for line in file:
            if line.startswith("AC"):
                acc = line.split()[1].split(';')[0]
            if line.startswith("DT") and "UniProtKB/Swiss-Prot" in line:
                sprot_ids.add(acc)

    return list(sprot_ids)


def map_and_extract_upkb_annots(
        parc2prot_map,
        uprot_annot_file,
        functions,
        out_file=None,
        no_return=True,
        uprot_idx=None,
        uparc_ids=None,
        sp_prefix="sp_",
        max_threads=2
):
    """ Specifically map Uniparc to UniprotKB annotations and extract. Swissprot sequences are additionally tagged. """

    if no_return and not out_file:
        raise RuntimeError("Merged annotations must either be returned as a dictionary, written to a file, or both.")

    # Get full Uniprot annot dictionary
    up_annots = map_and_extract_annots(
        parc2prot_map,
        uprot_annot_file,
        functions,
        uprot_idx=uprot_idx,
        uparc_ids=uparc_ids,
        record_type='uniprot',
        max_threads=max_threads
    )

    # Get Swissprot IDs, map to Uniparc IDs and extract annotations
    sp_ids = extract_sprot_ids(uprot_annot_file)
    prot2parc_map = {prot_id.split(',')[0] : parc_id
                     for parc_id, prot_id in file_util.tab_del_file_to_dict(parc2prot_map).items()
                     if prot_id}
    sp_ids_mapped = [prot2parc_map[sp_id] for sp_id in sp_ids]
    sp_annots = map_and_extract_annots(
        parc2prot_map,
        uprot_annot_file,
        functions,
        uprot_idx=uprot_idx,
        uparc_ids=sp_ids_mapped,
        record_type='uniprot',
        field_prefix=sp_prefix,
        max_threads=max_threads
    )

    merge_annots(
        annot_dicts=[up_annots, sp_annots],
        out_file=out_file,
        no_return=no_return
    )


def create_id_db(mapping_tsv, db_file):
    """ Creates an sqlite3 database of IDs from the second col of a tsv ID mapping file.
     NOTE: this is currently specific for the process of fast lookup during priority set generation for clustering. The
     purpose is simply fast lookup for the existence of at least one instance of the Uniparc ID mapped to the other
     database in the mapping file (e.g. SwissProt). """

    mapping_tsv = Path(mapping_tsv)
    db_file = Path(db_file)
    if db_file.exists():
        db_file.unlink()

    # Set up DB
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

    cur.execute("PRAGMA journal_mode=WAL;")
    cur.execute("PRAGMA synchronous=OFF;")
    cur.execute("PRAGMA temp_store=MEMORY;")

    # Just need a single table with single col
    cur.execute("DROP TABLE IF EXISTS mapped_ids")
    cur.execute("""
            CREATE TABLE mapped_ids (
                seq_id TEXT PRIMARY KEY
            )
        """)

    batch = []
    batch_size = 100000

    with open(mapping_tsv, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:  # Take just the second ID from each row
            seq_id = row[1].strip()
            batch.append((seq_id,))
            if len(batch) >= batch_size:
                # Just need the set of IDs - don't care about number of mappings
                cur.executemany(
                    "INSERT OR IGNORE INTO mapped_ids (seq_id) VALUES (?)",
                    batch
                )
                conn.commit()
                batch.clear()
        if batch:
            cur.executemany(
                "INSERT OR IGNORE INTO mapped_ids (seq_id) VALUES (?)",
                batch
            )
            conn.commit()

    conn.close()


def query_id_db(db_file, query_ids):
    """ Queries a list of IDs against an existing single table, single column sqlite DB and returns the subset of IDs
    present in the DB. """

    db_file = Path(db_file)
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

    # Get table name and ensure there's only one
    tables = cur.execute("""
            SELECT name
            FROM sqlite_master
            WHERE type='table'
              AND name NOT LIKE 'sqlite_%'
        """).fetchall()

    if len(tables) != 1:
        conn.close()
        raise ValueError("ID database should contain 1 table only.")

    table_name = tables[0][0]

    # Get col name and ensure there's only one
    col_info = cur.execute(f"PRAGMA table_info({table_name})").fetchall()
    if len(col_info) != 1:
        raise ValueError("ID database table should contain 1 column only.")

    col_name = col_info[0][1]

    # Create temp table containing query IDs
    cur.execute("CREATE TEMP TABLE query_ids (id TEXT PRIMARY KEY)")
    cur.executemany(
        "INSERT OR IGNORE INTO query_ids (id) VALUES (?)",
        ((x,) for x in query_ids)
    )

    # Get intersection of query IDs and original DB
    found = cur.execute(f"""
            SELECT q.id
            FROM query_ids q
            JOIN {table_name} m
              ON q.id = m.{col_name}
        """).fetchall()

    conn.close()

    # Return list of unique IDs in intersection
    return [row[0] for row in found]


def map_up_priority(
        uprot_annot_file=None,
        prot2parc_map=None,
        full_sprot_map_file=None,
        full_uprot_map_file=None,
        include=None
):
    """ Return Uniparc IDs mapping to Swissprot and all UniprotKB sequences as [[mapped_sprot], [mapped_uprot]]. If only
     a specific subset of Uniparc IDs are intended for inclusion, these can be specified. If include is None, all IDs
     present in the mapping file will be used. There is no problem having sequences in the include list which do not
     map to Uniprot sequences. Now updated to allow providing full SwissProt --> Uparc and UniProtKB --> Uparc mapping
     files """

    if (full_sprot_map_file and full_uprot_map_file and include):

        sprot_db_file = f"{full_sprot_map_file.rsplit('.', maxsplit=1)[0]}" \
                        f".sqlite"
        sprot_db_file = Path(sprot_db_file)

        uprot_db_file = f"{full_uprot_map_file.split('.', maxsplit=1)[0]}" \
                        f".sqlite"
        uprot_db_file = Path(uprot_db_file)

        # Create DB files if not already present
        if not sprot_db_file.exists():
            create_id_db(full_sprot_map_file, sprot_db_file)
        if not uprot_db_file.exists():
            create_id_db(full_uprot_map_file, uprot_db_file)

        # Get list of Uniparc query IDs mapped at least once to each DB
        sprot_mapped = query_id_db(sprot_db_file, include)
        uprot_mapped = query_id_db(uprot_db_file, include)

    elif (uprot_annot_file and prot2parc_map):

        prot2parc = file_util.tab_del_file_to_dict(prot2parc_map)

        if not include:  # If no subset is specified
            include = list(prot2parc.values())

        # Extract Swispprot IDs and map to Uniparc IDs
        sprot_ids = extract_sprot_ids(uprot_annot_file)
        sprot_mapped = [prot2parc[seq_id] for seq_id in sprot_ids if prot2parc[seq_id] in include]

        # Get all uniprot IDs
        uprot_mapped = [seq_id for seq_id in list(prot2parc.values()) if seq_id in include]

    else:
        raise RuntimeError(f"Either full Uniparc mapping files or dataset specific UniProt raw annotation + prot2parc "
                           f"mapping files must be provided.")

    return [sprot_mapped, uprot_mapped]


def extract_tax_lineage(tax_annots, tax_idx=None, uparc_ids=None, out_file=None, no_return=True, annot_merge=False):
    """ Extract all available lineage information for all or selected Uniparc entries.
    NOTE: """

    if no_return and not out_file:
        raise RuntimeError("If no return is specified, an annotation file name must be provided.")

    if not tax_idx:
        tax_idx = tax_annots.split('.')[0]+".idx"
        if tax_idx not in os.listdir():
            raise RuntimeError("Either a taxonomy annotation index file must be specified, or else an index file following "
                               "naming conventions of this module must be present in the directory.")

    if not uparc_ids:
        uparc_ids = [acc for acc in file_util.tab_del_file_to_dict(tax_idx).keys()]

    idx_dict = {
        acc : idx
        for acc, idx in file_util.tab_del_file_to_dict(tax_idx).items()
        if acc in uparc_ids
    }

    with open(tax_annots) as annot_file:

        tax_annot_dict = {}

        # TODO: Fix - tax_idx is just the name of the

        for acc, idx in idx_dict.items():

            tax_annot_dict[acc] = {}
            annot_file.seek(int(idx))
            xml_str = ""

            while True:
                line = annot_file.readline()
                xml_str += line
                if "</taxon>" in line:
                    break

            tax_dict = xmltodict.parse(xml_str)

            # Get lowest designation
            lowest_rank = tax_dict["taxon"]["@rank"]
            lowest_tax_id = tax_dict["taxon"]["@taxId"]
            lowest_tax = tax_dict["taxon"]["@scientificName"]
            low_desig_annot = f"{lowest_rank}: {lowest_tax}"  # Special annotation showing rank and name for lowest designation
            tax_annot_dict[acc]["tax_lowest"] = low_desig_annot
            # if tax_dict["taxon"]["@hidden"] == "false":
            #     tax_annot_dict[acc][f"tax_{lowest_rank}"] = lowest_tax
            #     tax_annot_dict[acc][f"tax_{lowest_rank}_id"] = lowest_tax_id
            # TODO: Fixing this
            tax_annot_dict[acc][f"tax_{lowest_rank}"] = lowest_tax
            tax_annot_dict[acc][f"tax_{lowest_rank}_id"] = lowest_tax_id

            # Get full lineage
            # NOTE: Currently collects only taxonomic designations for non-hidden ranks as well as any "clade"
            # designations (which are concatenated in a single string)
            for tax_level in tax_dict["taxon"]["lineage"]["taxon"]:

                if "@rank" in tax_level.keys() and tax_level["@hidden"] == "false":  # Include this annotation
                    tax_annot_dict[acc][f"tax_{tax_level['@rank']}"] = tax_level["@scientificName"]
                    tax_annot_dict[acc][f"tax_{tax_level['@rank']}_id"] = tax_level["@taxId"]

                elif "@rank" in tax_level.keys() and tax_level["@rank"] == "clade":
                    try:
                        tax_annot_dict[acc]["tax_clade"] += f";{tax_level['@scientificName']}"
                        tax_annot_dict[acc]["tax_clade_id"] += f";{tax_level['@taxId']}"
                    except KeyError:
                        tax_annot_dict[acc]["tax_clade"] = f"{tax_level['@scientificName']}"
                        tax_annot_dict[acc]["tax_clade_id"] = f"{tax_level['@taxId']}"

    if out_file:
        if not annot_merge:
            if out_file in os.listdir():
                raise RuntimeError(f"{out_file} already exists. To merge taxonomy annotations, annot_merge=True must "
                                   f"be specified.")
            else:
                annot_file_from_dict(out_file, tax_annot_dict)
        else:
            merge_annots(annot_dicts=tax_annot_dict, annot_files=out_file, out_file=out_file)

    if not no_return:
        return tax_annot_dict





def compound_annot(current_dict, label, function):
    """ From a current annotation dictionary, add an annotation which are dependent on one or more existing annotations.
    The annotation should be supplied as a lambda function, with references to existing annotations as:
    annots[seq][existing_label], where 'annots' and 'seq' will be defined variables within this function, and
    existing_label is the specific annotation being referenced. """

    annots = copy.deepcopy(current_dict)

    for seq in annots.keys():
        try:
            annots[seq][label] = function(annots, seq)
        except:   # If function fails, annotate as None
            annots[seq][label] = None

    annots = clean_annot_dict(annots)

    return annots



############################## PRE-DEFINED ANNOTATION VALUES, REFERENCES, FUNCTIONS ####################################

CLADE_DEFINE = [
        # Hard-coded (for now) clade definitions for superfamily tree
        # [name, in, out]

  ['RNaseZ', ['P39300', 'Q9WZW8'], 'Q9SID3'],
  ['RNaseJ', ['Q82ZZ3', 'Q72JJ7'], 'Q9SID3'],
  ['CPSF', ['Q9UKF6', '_NZ_APCS01000105.1_4'], 'Q9SID3'],
  ['DNA-repair', ['_ACE_DIBNCCHN_1_3532', 'Q96SD1'], 'Q9SID3'],
  ['B1/B2', ['P26918', 'P04190'], 'Q9SID3'],
  ['Alkyl-sulfatase', ['F8KAY7', 'Q9I5I9'], 'Q9SID3'],
  ['GLOXII', ['Q8ZRM2', 'Q9SID3'], 'Q72JJ7'],
  ['PSDO', ['Q9C8L4', 'C8WS08'], 'Q9SID3'],
  ['B3', ['Q89GW5', 'B5DCA0'], 'Q9SID3'],
  ['Lactonase', ['Q5W503', 'Q9X207'], 'Q9SID3'],
  ['FDP', ['Q9F0J6', 'P39695'], 'Q9SID3']
]


# class SelfRefDict:
#     """ Dictionary allowing self-referencing of other key:value pairs from the object. Used specifically here for
#      lambda annotation functions which may be inter-dependent. """
#
#     def __init__(self, dict:dict):
#         self.dict = dict
#
#     def __iter__(self):
#         return self.dict.__iter__()
#
#     def __getitem__(self, key):
#         return self.dict.__getitem__()
#
#     def internal_ref(self, label):
#         """ Refer to an internally defined annotation function within this object. """
#         return self[label]


############################### HELPER FUNCTIONS FOR UNIPROT ANNOTATION EXTRACTION #####################################

# Functions should take as input new-line-split Uniprot annotation record for a single ID

def get_up_seq(record_lines):
    """ Extract and return the sequence of a UnirpotKB record. """

    seq = ""
    in_seq = False
    for line in record_lines:
        if line.startswith("//"):  # End of record (and sequence)
            return seq
        elif line.startswith("SQ"):
            in_seq = True
        elif in_seq:  # Actual sequences starts on the line after the SQ flag
            seq += "".join(line.strip().split())

    if not in_seq:
        raise RuntimeError("No sequence found in provided annotation record.")
    else:
        raise RuntimeError("No end of record flag found in provided annotation record.")

########################################################################################################################

############################### HELPER FUNCTIONS FOR UNIPARC ANNOTATION EXTRACTION #####################################

# Functions should take as input Uniparc annotations for a single ID in XML ordered dictionary format

def uparc_interpro_hits(xml_dict):
    """ Extract all Interpro entry hit information from a Uniparc annotation record (XML format). """

    try:
        ipro_hits_dict = xml_dict["signatureSequenceMatch"]  # All interpro hits
    except KeyError:
        ipro_hits_dict = []

    # Return a list of tuples containing following info for each entry: (IPR ID, IPR Name, member_db, (start,end) )
    ipro_hits_tuple = []
    for hit in ipro_hits_dict:

        # Uniparc XML has format discrepancy when multiple hits are present for same Interpro tag
        if not isinstance(hit['lcn'], list):
            hit['lcn'] = [hit['lcn']]

        ipro_hits_tuple.append((hit['ipr']['@id'], hit['ipr']['@name'], hit['@database'],
                                [(region['@start'], region['@end']) for region in hit['lcn']]))

    return ipro_hits_tuple


# Collection of pre-defined lambda functions for querying Uniprot annotations
# Lambda variable should be a list (new-line-split) of the Uniprot annotation output format
UP_ANNOTS = {

    "sequence" : lambda x: get_up_seq(x),

    "accession" : lambda x: [line.split()[1].split(';')[0] for line in x if line.startswith('AC')][0],

    "length" : lambda x: [line.split(';')[1].split('AA')[0].strip() for line in x if line.startswith("ID")][0],

    "fragment" : lambda x: True if len([line for line in x if 'DE' in line and 'Fragment' in line]) > 0 else None,

    "organism" : lambda x: [line.split("OS")[1].strip() for line in x if line.startswith("OS")][0],

    "taxonomy" : lambda x: ' '.join([line.split("OC")[1].strip() for line in x if line.startswith("OC")]),

    "tax_id" : lambda x: [line.strip().split("NCBI_TaxID=")[1].split(";")[0] for line in x if line.startswith("OX")][0],

    "gene_name" : lambda x: [line.split('=')[1].split()[0] for line in x if line.startswith('GN')][0],

    "protein_name" : lambda x: [line.split("=")[1].split("{")[0].strip() for line in x if line.startswith("DE")][0],

    "ec_number" : lambda x: ';'.join(set([line.split("EC=")[1].split("{")[0].strip() for line in x if line.startswith("DE") and "EC=" in line])),

    "interpro_tag" : lambda x: ';'.join(sorted([line.split("; ")[1].strip() for line in x if 'InterPro;' in line])),

    "interpro_name" : lambda x: ';'.join(sorted([line.split("; ")[2].strip() for line in x if 'InterPro;' in line])),

    "cdd_tag" : lambda x: ';'.join(sorted([line.split("; ")[1].strip() for line in x if 'CDD;' in line])),

    "cdd_name" : lambda x: ';'.join(sorted([line.split("; ")[2].strip() for line in x if 'CDD;' in line])),

    "gene3d_tag" : lambda x: ';'.join(sorted([line.split("; ")[1].strip() for line in x if 'Gene3D;' in line])),

    "gene3d_name" : lambda x: ';'.join(sorted([line.split("; ")[2].strip() for line in x if 'Gene3D;' in line])),

    "panther_tag" : lambda x: ';'.join(sorted([line.split("; ")[1].strip() for line in x if 'PANTHER;' in line])),

    "panther_name" : lambda x: ';'.join(sorted([line.split("; ")[2].strip() for line in x if 'PANTHER;' in line])),

    "pfam_tag" : lambda x: ';'.join(sorted([line.split("; ")[1].strip() for line in x if 'Pfam;' in line])),

    "pfam_name" : lambda x: ';'.join(sorted([line.split("; ")[2].strip() for line in x if 'Pfam;' in line])),

    "pdb" : lambda x: ';'.join([line.split("; ")[1].strip() for line in x if 'PDB;' in line])

}

# Collection of pre-defined lambda functions for querying Uniparc annotations
# Lambda variable should be an XML ordered dictionary (as returned by uparc_xml_dict)
UPARC_ANNOTS = {

    "interpro_tag_parc" : lambda x: ';'.join(list(set([hit[0] for hit in uparc_interpro_hits(x)]))),

    "interpro_name_parc" : lambda x: ';'.join(list(set([hit[1] for hit in uparc_interpro_hits(x)])))

    # Uniprot ID mapping

}

COMPOUND_ANNOTS_UPROT = {

    # Interpro names excluding superfamily and generic MBL domain entries
    "interpro_specific_name" :
        lambda annots, seq: ';'.join([name for name in annots[seq]["interpro_name"].split(';') if name not in ["Metallo-B-lactamas.", "RibonucZ/Hydroxyglut_hydro."]])

}