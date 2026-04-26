import atexit
import gzip
import os
import random
import shutil
import subprocess
import threading
import multiprocessing as mp
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed, wait
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from . import hmm
from . import file_util
from . import annots

""" Functions and workflows for profile-based curation from sequence databases.  """


def fetch_uniparc_fasta(
        segments=None,
        path=None,
        unzip=True,
        max_threads=1,
        verbose=1
):
    """ Download one, multiple or all volumes of current Uniparc release in fasta (gzip) format. Desired volumes should
     be specified as integers. """

    if segments and not isinstance(segments, list):
        segments = [segments]
    elif not segments:
        # NOTE: This is hard-coded as number of segments from previous few releases (=200)
        segments = range(1, 201)

    def get_segment(segment, path=None, unzip=False, verbose=verbose):

        segment_name = f"uniparc_active_p{segment}.fasta"
        zipped_name = segment_name + ".gz"
        address = f"https://ftp.expasy.org/databases/uniprot/current_release/uniparc/fasta/active/{zipped_name}"
        
        if verbose == 1:
            print(f"Commencing download for segment {segment}.")

        wget_cmd = ["wget"]

        if verbose != 2:
            wget_cmd.append("--quiet")

        if path:
            wget_cmd.extend(["-O", path])
            seg_p = Path(path)/segment_name
            zipped_p = Path(path)/zipped_name
        else:
            seg_p = Path(segment_name)
            zipped_p = Path(zipped_name)

        wget_cmd.append(address)
        subprocess.run(wget_cmd, check=True)
        
        if unzip:
            with gzip.open(zipped_p, 'rb') as in_f: 
                with open(Path(seg_p), 'wb') as out_f: 
                    shutil.copyfileobj(in_f, out_f)

        if verbose == 1:
            print(f"Finished download for segment {segment}.")

    # Submit download jobs as threads
    with ThreadPoolExecutor(max_workers=max_threads) as executor:

        try:
            future_list = [executor.submit(get_segment, segment, path, unzip) for segment in segments]
            for future in as_completed(future_list):
                future.result()
        
        finally:
            executor.shutdown(wait=False, cancel_futures=True)


def profile_uniparc_search(
        profile,
        hit_file,
        threshold=0,
        best_domain=True,
        max_download_threads=1,
        max_hmm_threads=1,
        nice=None,
        heuristics_off=False,
        remove_segments=True
):
    """ Search Uniparc against a profile HMM extracting hits above a specified bit-score threshold to a single fasta
     file. Max_threads represents the TOTAL maximum number of threads to use, distributed equally over max_processes. """

    lock = threading.Lock()

    # Uniparc segments
    segments = range(1, 201)

    # Evenly divide available hmm threads among parallel segment searches
    threads_per_segment = max_hmm_threads // max_download_threads

    def profile_search_segment(segment):

        fetch_uniparc_fasta(segment, unzip=True)  # Download and unzip this segment

        # Profile search for this segment
        score_dict = hmm.hmm_search(
            profile,
            f"uniparc_active_p{segment}.fasta",
            results_file=f"temp_results_{segment}.txt",
            id_format=None,
            best_domain=best_domain,
            max_on=heuristics_off,
            max_threads=threads_per_segment,
            nice=nice
        )

        hits_above = []

        for seq, score in score_dict.items():   # Check which sequences score above threshold
            if score > threshold:
                hits_above.append(seq)

        # Extract hits above threshold from fasta
        file_util.raw_extract_fasta(
            f"uniparc_active_p{segment}.fasta",
            f"temp_{segment}.fasta",
            target_seqs=hits_above,
            id_format=None
        )

        with lock:
            os.system(f"cat temp_{segment}.fasta >> {hit_file}")

        # Remove temp files
        os.remove(f"temp_{segment}.fasta")
        os.remove(f"temp_results_{segment}.txt")
        if remove_segments:
            os.remove(f"uniparc_active_p{segment}.fasta")

    # Multiple processes
    executor = ThreadPoolExecutor(max_workers=max_download_threads)
    print("Thread pool created")    # TODO: For de-bugging - DELETE
    future_list = [executor.submit(profile_search_segment, segment) for segment in segments]
    atexit.register(executor.shutdown, wait=False)
    print(len(future_list))  # TODO: For de-bugging - DELETE
    print("Jobs submitted")  # TODO: For de-bugging - DELETE
    wait(future_list)


def profile_uniparc_search_new(
        profiles,
        hit_file,
        thresholds,
        dom_hit_file=None,
        uparc_file_path=None,
        best_chain=False,
        best_domain=True,
        max_download_threads=1,
        max_hmm_threads=1,
        nice=None,
        max_on=False,
        remove_segments=True,
        region_type="env",
        overlap_tol_prop=0.1
):
    """ Search Uniparc against a profile HMM extracting hits above a specified bit-score threshold to a single fasta
     file. Max_threads represents the TOTAL maximum number of threads to use, distributed equally over max_processes. 
     
     TODO: Implement Extraction of subsequences on basis of hit coordiantes
     TODO: (will do this when I re-format/standardise profile_db etc.) """

    if not isinstance(profiles, list):  # Single profile only
        profiles = [profiles]

    if not isinstance(thresholds, list):
        thresholds = [thresholds]

    if len(profiles) != len(thresholds):
        raise RuntimeError("The number of provided profiles must match number of provided thresholds.")

    if uparc_file_path:
        for file in [file for file in os.listdir(uparc_file_path) if "uniparc_active_p" in file]:
            subprocess.run(["ln", "-s", uparc_file_path.rstrip("/") + f"/{file}", "."])

    lock = mp.Lock()

    # Uniparc segments - TODO: hard-coded for # of segments (stable for several releases now)
    # segments = range(1, 201)
    segments = range(1,17)  # TODO: TESTING ON REDUCED UNIPARC DB - REVERT

    # Evenly divide available hmm threads among parallel segment searches
    threads_per_segment = max_hmm_threads // max_download_threads

    # def search_segment(segment):
    #
    #     if not (os.path.exists(f"uniparc_active_p{segment}.fasta") or os.path.islink(f"uniparc_active_p{segment}.fasta")):  # Check if segment already present, else download it
    #         fetch_uniparc_fasta(segment, unzip=True)  # Download and unzip this segment
    #
    #     if best_chain:
    #
    #         best_chain_scores = hmm.hmm_search_chain(profiles, f"uniparc_active_p{segment}.fasta",
    #                                                  results_file=f"temp_results_{segment}.txt", region_type=region_type,
    #                                                  overlap_tol_prop=overlap_tol_prop, max_on=max_on,
    #                                                  max_threads=threads_per_segment, nice=nice)
    #
    #         hits_above = set()
    #         for i in range(len(profiles)):
    #             for target, score in best_chain_scores[profiles[i].split(".hmm")[0]].items():
    #                 if score >= thresholds[i]:
    #                     hits_above.add(target)
    #
    #         parsers.extract_fasta(f"uniparc_active_p{segment}.fasta", f"temp_{segment}.fasta", list(hits_above))
    #
    #     else:
    #
    #         if len(profiles) == 1:  # Single profile - use hmm_search
    #
    #             # Profile search for this segment
    #             score_dict = hmm.hmm_search(profiles[0], f"uniparc_active_p{segment}.fasta",
    #                                         results_file=f"temp_results_{segment}.txt",
    #                                         id_format=None, best_domain=best_domain, max_on=max_on,
    #                                         max_threads=threads_per_segment, nice=nice)
    #
    #             hits_above = []
    #
    #             for seq, score in score_dict.items():   # Check which sequences score above threshold
    #                 if score > thresholds[0]:
    #                     hits_above.append(seq)
    #
    #             # Extract hits above threshold from fasta
    #             parsers.raw_extract_fasta(f"uniparc_active_p{segment}.fasta", f"temp_{segment}.fasta", target_seqs=hits_above,
    #                                     id_format=None)
    #
    #         else:  # Multiple profiles - use profile_db
    #
    #             hmm.profile_db(profiles, thresholds, f"uniparc_active_p{segment}.fasta", hits_file=f"temp_{segment}.fasta",
    #                         name=f"all_{segment}", results_file=f"temp_results_{segment}.txt",
    #                         best_domain=best_domain, max_threads=threads_per_segment, nice=nice)
    #
    #             os.remove(f"all_{segment}.hmm")
    #
    #
    #     with lock:
    #         # FOR TESTING - DELETE
    #         print(f"Writing hits for segment {segment}")
    #         os.system(f"cat temp_{segment}.fasta >> {hit_file}")
    #
    #     # Remove temp files
    #     os.remove(f"temp_{segment}.fasta")
    #     os.remove(f"temp_results_{segment}.txt")
    #     if remove_segments:
    #         os.remove(f"uniparc_active_p{segment}.fasta")


    # # TODO: FOR TESTING - DELETE
    # profile_search_segment(single_segment)


    # TODO: UNCOMMENT BELOW 8 LINES - FOR TESTING

    # Multiple processes
    ctx = mp.get_context("spawn")
    executor = ProcessPoolExecutor(max_workers=max_download_threads, mp_context=ctx)
    print("Process pool created")    # TODO: For de-bugging - DELETE
    future_list = [executor.submit(search_segment, segment, profiles, thresholds, best_chain, best_domain,
                                   region_type, max_on, overlap_tol_prop, threads_per_segment, nice, remove_segments)
                   for segment in segments]

    # atexit.register(executor.shutdown, wait=False)

    print(f"{len(future_list)} jobs submitted.")  # TODO: For de-bugging - DELETE
    try:
        for fut in as_completed(future_list):
            fut.result()
    except Exception:
        for f in future_list:
            f.cancel()
        # executor.shutdown(wait=False)
        raise
    except KeyboardInterrupt:
        for f in future_list:
            f.cancel()
        # executor.shutdown(wait=False)
        raise
    finally:
        executor.shutdown(wait=False)

    print("Concatenating full segment files")
    if hit_file in os.listdir():
        os.remove(hit_file)
    with open(hit_file, "a") as full_file:
        for segment in segments:
            with open(f"temp_{segment}.fasta") as seg_file:
                shutil.copyfileobj(seg_file, full_file)
            os.remove(f"temp_{segment}.fasta")

    if dom_hit_file and best_chain:  # TODO: Coords should soon be implemented for non-chain searches too
        print("Concatenating domain segment files")
        if dom_hit_file in os.listdir():
            os.remove(dom_hit_file)
        with open(dom_hit_file, "a") as full_dom_file:
            for segment in segments:
                with open(f"temp_{segment}_dom.fasta") as dom_seg_file:
                    shutil.copyfileobj(dom_seg_file, full_dom_file)
                os.remove(f"temp_{segment}_dom.fasta")


def search_segment(
        segment,
        profiles,
        thresholds,
        best_chain,
        best_domain,
        region_type,
        max_on,
        overlap_tol_prop,
        threads_per_segment,
        nice,
        remove_segments
):

    if not (os.path.exists(f"uniparc_active_p{segment}.fasta")
            or os.path.islink(f"uniparc_active_p{segment}.fasta")):  # Check if segment already present, else download it
        fetch_uniparc_fasta(segment, unzip=True)  # Download and unzip this segment

    if best_chain:

        best_chain_scores = hmm.hmm_search_chain(
            profiles,
            f"uniparc_active_p{segment}.fasta",
            results_file=f"temp_results_{segment}.txt",
            region_type=region_type,
            overlap_tol_prop=overlap_tol_prop,
            max_on=max_on,
            max_threads=threads_per_segment,
            nice=nice)

        hits_above = []
        hits_above_coords = []  # Coordinates of chain
        for i in range(len(profiles)):
            for target, best_chain in best_chain_scores[profiles[i].split(".hmm")[0]].items():
                if best_chain[0] >= thresholds[i]:
                    hits_above.append(target)
                    hits_above_coords.append(best_chain[1])  # For extracting domain subseqs

        # TODO: Need to determine chain with longest span per sequence - currently domain coords will correspond to first profile only

        # Generate file of full hit seqs and domain subseqs from hits
        file_util.extract_fasta(
            f"uniparc_active_p{segment}.fasta", f"temp_{segment}.fasta", hits_above
        )

        file_util.extract_fasta(
            f"uniparc_active_p{segment}.fasta",
            f"temp_{segment}.fasta",
            hits_above)

        file_util.extract_fasta(f"uniparc_active_p{segment}.fasta",
                                f"temp_{segment}_dom.fasta",
                                hits_above,
                                coords=hits_above_coords)

    else:

        if len(profiles) == 1:  # Single profile - use hmm_search

            # Profile search for this segment
            score_dict = hmm.hmm_search(
                profiles[0],
                f"uniparc_active_p{segment}.fasta",
                results_file=f"temp_results_{segment}.txt",
                id_format=None,
                best_domain=best_domain,
                max_on=max_on,
                max_threads=threads_per_segment,
                nice=nice
            )

            hits_above = []

            for seq, score in score_dict.items():  # Check which sequences score above threshold
                if score > thresholds[0]:
                    hits_above.append(seq)

            # Extract hits above threshold from fasta
            file_util.raw_extract_fasta(
                f"uniparc_active_p{segment}.fasta",
                f"temp_{segment}.fasta",
                target_seqs=hits_above,
                id_format=None)

        else:  # Multiple profiles - use profile_db

            hmm.profile_db(
                profiles,
                thresholds,
                f"uniparc_active_p{segment}.fasta",
                hits_file=f"temp_{segment}.fasta",
                name=f"all_{segment}",
                results_file=f"temp_results_{segment}.txt",
                max_on=max_on,
                best_domain=best_domain,
                max_threads=threads_per_segment,
                nice=nice)

            os.remove(f"all_{segment}.hmm")

    print(f"Segment {segment} scan completed")

    # with lock:  # TODO: Now performing this step in parent process - remove here once tested
    #     # FOR TESTING - DELETE
    #     print(f"Writing hits for segment {segment}")
    #     os.system(f"cat temp_{segment}.fasta >> {hit_file}")

    # Remove temp files
    # os.remove(f"temp_{segment}.fasta")  # TODO: Now performing this step in parent process - remove here once tested
    os.remove(f"temp_results_{segment}.txt")
    if remove_segments:
        os.remove(f"uniparc_active_p{segment}.fasta")


def profile_search_segment_refseq(
        profiles,
        hit_file,
        thresholds,
        db_dir=None,
        segment_files=None,
        best_domain=True,
        max_download_threads=1,
        max_hmm_threads=1,
        nice=None
):
    """ Search a segmented database against a profile HMM extracting hits above a specified bit-score threshold to a single fasta
     file. Max_threads represents the TOTAL maximum number of threads to use, distributed equally over max_processes. """

    # NOTE: Retrieval of segments is not implemented here - this function is generic for any segmented database
    #
    if not isinstance(profiles, list):  # Single profile only
        profiles = [profiles]

    if not isinstance(thresholds, list):
        thresholds = [thresholds]

    if len(profiles) != len(thresholds):
        raise RuntimeError("The number of provided profiles must match number of provided thresholds.")

    if not segment_files:
        if db_dir:  # Assume we want to run over all segment files - get full list
            segment_files = [f"{db_dir}/{file}" for file in os.listdir(db_dir) if ".fa" in file]
        else:
            raise RuntimeError("Either a directory containing all database segments or list of segment file paths must "
                               "be provided.")

    lock = threading.Lock()

    # Evenly divide available hmm threads among parallel segment searches
    threads_per_segment = max_hmm_threads // max_download_threads

    def search_segment(segment):

        seg_file_name = segment.split("/")[-1]
        subprocess.run(["ln", "-s", segment, "temp_" + seg_file_name])

        hmm.profile_db(
            profiles,
            thresholds,
            "temp_"+seg_file_name,
            hits_file=f"temp_hits_{seg_file_name}.fasta",
            name=f"all_{seg_file_name}",
            results_file=f"temp_results_{seg_file_name}.txt",
            best_domain=best_domain,
            max_threads=threads_per_segment,
            nice=nice)

        if len(profiles) > 1:  # A collated .hmm file will have been created by profile_db
            os.remove(f"all_{seg_file_name}.hmm")

        with lock:
            # FOR TESTING - DELETE
            print(f"Writing hits for segment {seg_file_name}")
            os.system(f"cat temp_hits_{seg_file_name}.fasta >> {hit_file}")

        # Remove temp files
        os.remove(f"temp_{seg_file_name}")
        os.remove(f"temp_hits_{seg_file_name}.fasta")
        os.remove(f"temp_results_{seg_file_name}.txt")

    # Multiple threads
    executor = ThreadPoolExecutor(max_workers=max_download_threads)
    print("Thread pool created")    # TODO: For de-bugging - DELETE
    future_list = [executor.submit(search_segment, segment) for segment in segment_files]
    atexit.register(executor.shutdown, wait=False)
    print(len(future_list))  # TODO: For de-bugging - DELETE
    print("Jobs submitted")  # TODO: For de-bugging - DELETE
    wait(future_list)


def profile_refseq_search(
        profiles,
        thresholds,
        refseq_dir,
        hit_file,
        best_domain=True,
        max_download_threads=1,
        max_hmm_threads=1,
        nice=None
):

    """
    # TODO: Will need to be changed to use new generic segmented database search funciton

    NOTE: Currently assumes that RefSeq database is downloaded and available in current directory. """

    # TODO: I think everything dowsntream should handle gzipped files now but leaving here for now just in case
    # # Segments will need to be unzipped for various steps
    # to_unzip = [f"{refseq_dir}/{file}" for file in os.listdir(refseq_dir) if (".gz" in file)]
    # with futures.ThreadPoolExecutor(max_workers=max_hmm_threads) as executor:
    #     future_list = [executor.submit(lambda seg: subprocess.run(["gzip", "-d", seg]), segment)
    #                    for segment in to_unzip]
    #     atexit.register(executor.shutdown, wait=False)

    profile_search_segment_refseq(profiles, hit_file, thresholds, refseq_dir, best_domain=best_domain,
                           max_download_threads=max_download_threads, max_hmm_threads=max_hmm_threads, nice=nice)

    # # Re-gzip any files that were unzipped
    # for file in [file.split(".gz")[0] for file in to_unzip]:
    #     subprocess.run(["gzip", file])


def pss_init(
        profiles,
        thresholds,
        hmm_mode,
        region_type,
        chain_overlap_tol,
        max_on,
        remove_seg_files,
        idx_suffix,
        cpu_per_search,
        nice
):
    """ Process initialiser for profile_search_segment. """

    global _profiles
    _profiles = profiles

    global _thresholds
    _thresholds = thresholds

    global _hmm_mode
    _hmm_mode = hmm_mode

    global _region_type
    _region_type = region_type

    global _chain_overlap_tol
    _chain_overlap_tol = chain_overlap_tol

    global _max_on
    _max_on = max_on

    global _remove_seg_files
    _remove_seg_files = remove_seg_files

    global _idx_suffix
    _idx_suffix = idx_suffix

    global _cpu_per_search
    _cpu_per_search = cpu_per_search

    global _nice
    _nice = nice


def profile_search_segment(
        segment,
        idx_file,
        seg_hits_file,
        seg_hits_dom_file
):
    """ Per-process searching of a single segment using profile searching with the specified hmm mode. """

    hmm_results = hmm.hmm_search_custom(
        _profiles,
        segment,
        mode=_hmm_mode,
        region_type=_region_type,
        chain_overlap_tol=_chain_overlap_tol,
        max_on=_max_on,
        cpu_per_search=_cpu_per_search,
        nice=_nice
    )

    hits = set()
    coords = {}

    # Apply thresholds
    for i in range(len(_profiles)):
        profile_name = _profiles[i].split(".hmm")[0]
        profile_threshold = _thresholds[i]
        for seq, result in hmm_results[profile_name].items():
            if result[0] >= profile_threshold:
                hits.add(seq)
                try:
                    coords[seq].append(result[1])
                except KeyError:
                    coords[seq] = [result[1]]

    # Extract domain based on lower median boundary over all profile by hit length
    hits = list(hits)
    median_bounds = []
    for hit in hits:
        
        # Sort low to high by total hit length
        hit_bounds = coords[seq]
        hit_bounds.sort(key=lambda x: x[1]-x[0])
        n = len(hit_bounds)
        median_bounds.append(hit_bounds[(n-1)//2])

    # Infer idx file if suffix provided
    if _idx_suffix and not idx_file:
        idx_file = segment.with_suffix("." + _idx_suffix.lstrip('.'))

    # Extract hits to file
    file_util.extract_fasta(
        segment,
        seg_hits_file,
        target_seqs=hits,
        idx_file=idx_file
    )

    if seg_hits_dom_file:
        file_util.extract_fasta(
            segment,
            seg_hits_dom_file,
            target_seqs=hits,
            idx_file=idx_file,
            coords=median_bounds
        )


def profile_search_segmented_db(
        profiles,
        thresholds,
        seg_files,
        hit_file,
        seg_idx_files=None,
        idx_suffix=None,
        dom_hit_file=None,
        hmm_mode="dom",
        region_type="env",
        chain_overlap_tol=0.1,
        max_on=False,
        remove_seg_files=False,
        cpu_per_search=1,
        max_search_processes=6,
        nice=None
):
    """ Generic function for multi-process-enabled pHMM searching of any segmented database (fasta format only)

     TODO: This should replace the Uniparc and RefSeq specific implementations above
     """

    if hit_file in os.listdir():
        os.remove(hit_file)
    if dom_hit_file and dom_hit_file in os.listdir():
        os.remove(dom_hit_file)

    # If idx files are provide explicitly
    if seg_idx_files:

        if len(seg_files) != len(seg_idx_files):
            raise ValueError("An index file must be provided for each segment "
                            "if provided explicitly.")
        
        idx_map = {seg_files[i] : seg_idx_files[i]
                   for i in range(len(seg_files))}
        
    else:
        idx_map = {seg_files[i] : None
                   for i in range(len(seg_files))}

    temp_tag = random.randint(int(1e8), int(1e9)-1)

    seg_hit_files = {
        file : f"{temp_tag}_{str(Path(file.name).with_suffix('.fa'))}"
        for file in seg_files
    }

    dom_seg_hit_files = {
        # TODO: Need to update this when moving to per-iteration output directories
        file :
            seg_hit_files[file].split(".")[0]+"_dom.fa"
            if dom_hit_file else None
        for file in seg_files
    }

    with ProcessPoolExecutor(
            max_workers=max_search_processes,
            initializer=pss_init,
            initargs=(
                    profiles,
                    thresholds,
                    hmm_mode,
                    region_type,
                    chain_overlap_tol,
                    max_on,
                    remove_seg_files,
                    idx_suffix,
                    cpu_per_search,
                    nice
            )
    ) as executor:

        # Submit search jobs
        future_seg_map = {executor.submit(profile_search_segment,
                                          segment,
                                          idx_map[segment],
                                          seg_hit_files[segment],
                                          dom_seg_hit_files[segment])
                          : segment
                          for segment in seg_files}

        # Write segment hits to file/s once search complete
        for future in as_completed(future_seg_map):

            try:
                future.result()
            except Exception:
                for f in future_seg_map:
                    f.cancel()
                raise

            segment = future_seg_map[future]
            seg_hit_file = seg_hit_files[segment]

            with open(hit_file, "a") as hit_f:
                with open(seg_hit_file) as seg_hit_f:
                    shutil.copyfileobj(seg_hit_f, hit_f)

            os.remove(seg_hit_file)

            if dom_hit_file:
                dom_seg_hit_file = dom_seg_hit_files[segment]
                with open(dom_hit_file, "a") as dom_hit_f:
                    with open(dom_seg_hit_file) as dom_seg_hit_f:
                        shutil.copyfileobj(dom_seg_hit_f, dom_hit_f)

                os.remove(dom_seg_hit_file)


def fetch_ena_proteins(seq_ids, out_file=None, no_return=True):

    if no_return and not out_file:
        raise ValueError(
            "If no_return is specified, an output fasta " "file must be provided"
        )
    
    records = []

    # Batch into chunks of 200
    batches = []
    for i in range(0, len(seq_ids), 200):
        batches.append(seq_ids[i : i + 200])

    # Fetch batches and process
    for batch in batches:

        result = annots.fetch_batch_raw(
            batch,
            dbName="ena_coding",
            format="default")[0]
        
        entries = [entry for entry in result.split("//\n")
                   if "ID" in entry]
        for entry in entries:
            prot_id = entry.split('/protein_id="')[1].split('"')[0]
            seq_str = ""
            found_prot = False
            for line in entry.splitlines():
                if "/translation=" in line:
                    found_prot = True
                    seq_str += line.strip().split('/translation="')[1]
                elif found_prot:
                    seq_str += line.strip().split()[1].strip('"')
                    if '"' in line:
                        break
            if prot_id not in [seq.id for seq in records]:
                records.append(SeqRecord(
                                    seq=Seq(seq_str),
                                    id=prot_id,
                                    description=""))
            
    if out_file:
        SeqIO.write(records, out_file, "fasta")
    if not no_return:
        return records

