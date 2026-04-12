import atexit
import os
import random
import shutil
import subprocess
import threading
import multiprocessing as mp
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed, wait
from pathlib import Path

from . import hmm
from . import file_util

""" Functions and workflows for profile-based curation from sequence databases.  """


def fetch_uniparc_fasta(
        segments=None,
        unzip=False,
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

    def get_segment(segment, unzip=False, verbose=verbose):
        # address = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniparc/fasta/active/uniparc_active_p{segment}.fasta.gz"
        address = f"https://ftp.expasy.org/databases/uniprot/current_release/uniparc/fasta/active/uniparc_active_p{segment}.fasta.gz"
        if verbose == 1:
            print(f"Commencing download for segment {segment}.")
        os.system(f"wget {'--quiet ' if verbose != 2 else ''}{address}")
        if unzip:
            os.system(f"gzip -d uniparc_active_p{segment}.fasta.gz")
        if verbose == 1:
            print(f"Finished download for segment {segment}.")

    # Multi-threading
    executor = ThreadPoolExecutor(max_workers=max_threads)
    future_list = [executor.submit(get_segment, segment, unzip) for segment in segments]
    wait(future_list)


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
        subprocess.run(["ln", "-s", segment, "temp_"+seg_file_name])

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
        seg_hits_file,
        seg_hits_dom_file
):
    """ Per-process searching of a single segment using profile searching with the specified hmm mode. """

    # TODO: Soon to implement a single call point for phmm searching which returns (score, coords) for any mode

    # TODO: Currently just implementing chain searching - only one that currently meets return contract

    testing_tag = random.randint(int(1e8), int(1e9-1)) # TODO: For optimisation testing - DELETE
    test_log_file = f"{testing_tag}_seg_{segment}.log"
    with open(test_log_file, 'w') as log_f:
        pass

    start = time.time()

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

    end = time.time()

    with open(test_log_file, 'a') as log_f:
        log_f.write(f"HMM searching of segment completed in {round(end-start)} seconds.\n\n")

    hits = set()
    coords = {}

    start = time.time()

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

    end = time.time()

    with open(test_log_file, 'a') as log_f:
        log_f.write(f"Application of thresholds to hits completed in {round(end-start)} seconds.\n\n")

    start = time.time()

    # Extract earliest profile start and latest profile end
    hits = list(hits)
    widest_bounds = [  # Overall region to extract for each seq
        (min(coords[seq], key=lambda x: x[0])[0],
         max(coords[seq], key=lambda x: x[1])[1])
        for seq in hits
    ]

    end = time.time()

    if _idx_suffix:
        idx_file = segment.split('.')[0] + '.' + _idx_suffix.lstrip('.')
    else:
        idx_file = None

    with open(test_log_file, 'a') as log_f:
        log_f.write(f"Determination of extraction endpoints completed in {round(end-start)} seconds.\n\n")

    start = time.time()

    # Extract hits to file
    file_util.extract_fasta(
        segment,
        seg_hits_file,
        target_seqs=hits,
        idx_file=idx_file
    )

    end = time.time()

    with open(test_log_file, 'a') as log_f:
        log_f.write(f"Extraction of full hits completed in {round(end-start)} seconds.\n\n")

    start = time.time()

    if seg_hits_dom_file:
        file_util.extract_fasta(
            segment,
            seg_hits_dom_file,
            target_seqs=hits,
            idx_file=idx_file,
            coords=widest_bounds
        )

    end = time.time()

    with open(test_log_file, 'a') as log_f:
        log_f.write(f"Extraction of domain hits completed in {round(end-start)} seconds.\n\n")

    Path(test_log_file).unlink()


def profile_search_segmented_db(
        profiles,
        thresholds,
        seg_files,
        hit_file,
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
        os.remove(dom_hit_file)
    if dom_hit_file and dom_hit_file in os.listdir():
        os.remove(dom_hit_file)

    temp_tag = random.randint(int(1e8), int(1e9)-1)

    seg_hit_files = {
        file : f"{temp_tag}_{file.split('.')[0]}.fa"
        for file in seg_files
    }

    dom_seg_hit_files = {
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
                                          seg_hit_files[segment],
                                          dom_seg_hit_files[segment])
                          : segment
                          for segment in seg_files}

        # Write segment hits to file/s once search complete
        for future in as_completed(future_seg_map):

            future.result()
            segment = future_seg_map[future]
            seg_hit_file = seg_hit_files[segment]

            with open(hit_file, 'a') as hit_f:
                with open(seg_hit_file) as seg_hit_f:
                    shutil.copyfileobj(seg_hit_f, hit_f)

            os.remove(seg_hit_file)

            if dom_hit_file:
                dom_seg_hit_file = dom_seg_hit_files[segment]
                with open(dom_hit_file, 'a') as dom_hit_f:
                    with open(dom_seg_hit_file) as dom_set_hit_f:
                        shutil.copyfileobj(dom_set_hit_f, dom_hit_f)

                os.remove(dom_seg_hit_file)

