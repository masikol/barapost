# *-* coding: utf-8 *-*
# In this module, kernel functions for prober are defined.

import os

from src.printlog import printl
from src.filesystem import remove_tmp_files

from src.prune_seqs import prune_seqs
from src.write_classification import write_classification
from src.prober_modules.networking import configure_request, send_request, wait_for_align, BlastError
from src.prober_modules.prober_spec import parse_align_results_xml, write_hits_to_download

from src.fasta import fasta_packets_from_str

def _split_and_resubmit(packet, packet_size, packet_mode, pack_to_send, seqs_processed,
        fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
        blast_algorithm, user_email, organisms,
        acc_dict, out_of_n):
    """
    :param packet: "packet" dictionary described in "barapost-prober.py" before the kernel loop:
    :type packet: dict;
    :param packet_size: size of the packet (see option '-c' for definition);
    :type packet_size: int;
    :param packet_mode: packet forming mode (see option '-c' for definition);
    :type packet_mode: int;
    :param pack_to_send: ordinal number of packet to send (it is list rather that in because it should be mutable);
    :type pack_to_send: list<int>;
    :param seqs_processed: nuber of sequnces processed (it is list rather that in because it should be mutable);
    :type seqs_processed: list<int>;
    :param fq_fa_path: path to current input file;
    :type fq_fa_path: str;
    :param tmp_fpath: path to current temporary file;
    :type tmp_fpath: str;
    :param taxonomy_path: path to taxonomt file;
    :type taxonomy_path: str;
    :param tsv_res_path: path to current classification file;
    :type tsv_res_path: str;
    :param acc_fpath: path to file 'hits_to_download.tsv';
    :type acc_fpath: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    :param blast_algorithm: BLAST algorithm to use (see option '-a');
    :type blast_algorithm: str;
    :param user_email: user email ot send with request;
    :type user_email: str;
    :param organisms: list of strings performing 'nt' database slices;
    :type organisms: list<str>;
    :param acc_dict: accession dictionary for writing to 'hits_to_download.tsv';
    :type acc_dict: dict<str: (str, int)>;
    :param out_of_n: dictionary for printing how many packets left;
    :type out_of_n: dict<str: str, str: int>;
    """

    # Number of sequnces in packet to be splitted:
    pack_len = len(packet["qual"])

    if pack_len > 1:
        # Split current packet into two (of equal number of sequences) and resubmit them one-by-one

        printl(logfile_path, "Splitting current packet into two and submitting each of them one-by-one.")

        # Update this dictionary to print how many packets left
        if not out_of_n["npacks"] is None:
            out_of_n["npacks"] += 1
            out_of_n["msg"] = " out of {}".format(out_of_n["npacks"])
        # end if

        # Calculate size of subpacket
        new_pack_size_0 = pack_len // 2
        if pack_len % 2 != 0:
            new_pack_size_0 += 1
        # end if

        # Split the packet
        for i, splitted_packet in enumerate(fasta_packets_from_str(packet["fasta"], new_pack_size_0)):

            # Inherit quality information from "ancestor" qual_dict
            for query_name in splitted_packet["qual"].keys():
                splitted_packet["qual"][query_name] = packet["qual"][query_name]
            # end for

            # Submit subpacket
            submit(splitted_packet, new_pack_size_0, 0, pack_to_send, seqs_processed,
                fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
                blast_algorithm, user_email, organisms,
                acc_dict, out_of_n)
        # end for
    else:
        # Prune the only sequence in packet and resend it

        printl(logfile_path, "Current packet contains only one sequence.")
        printl(logfile_path, "prober will prune this sequence twofold and resubmit it.")

        # Calculate new length for this sequence
        old_seq = map(str.strip, packet["fasta"].splitlines()[1:]) # generator of stripped sequence-sontaining lines
        old_len = len(''.join(old_seq)) # calculate length of old sequence
        new_len = old_len // 2
        if old_len % 2 != 0:
            new_len += 1
        # end if

        packet["fasta"] = prune_seqs(packet["fasta"], new_len)

        submit(packet, packet_size, packet_mode, pack_to_send, seqs_processed,
            fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
            blast_algorithm, user_email, organisms,
            acc_dict, out_of_n)
    # end if
# end def _split_and_resubmit


def _handle_result(align_xml_text, packet, taxonomy_path,
    tsv_res_path, acc_dict, acc_fpath, seqs_processed, pack_to_send,
    tmp_fpath, logfile_path):
    """
    :param align_xml_text: XML text with results of alignment;
    :type align_xml_text: str;
    :param packet: "packet" dictionary described in "barapost-prober.py" before the kernel loop:
    :type packet: dict;
    :param taxonomy_path: path to taxonomt file;
    :type taxonomy_path: str;
    :param tsv_res_path: path to current classification file;
    :type tsv_res_path: str;
    :param acc_dict: accession dictionary for writing to 'hits_to_download.tsv';
    :type acc_dict: dict<str: (str, int)>;
    :param acc_fpath: path to file 'hits_to_download.tsv';
    :type acc_fpath: str;
    :param seqs_processed: nuber of sequnces processed (it is list rather that in because it should be mutable);
    :type seqs_processed: list<int>;
    :param pack_to_send: ordinal number of packet to send (it is list rather that in because it should be mutable);
    :type pack_to_send: list<int>;
    :param tmp_fpath: path to current temporary file;
    :type tmp_fpath: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    """

    # Get result tsv lines
    result_tsv_lines = parse_align_results_xml(align_xml_text,
        packet["qual"], acc_dict, logfile_path, taxonomy_path)

    # Write classification to TSV file
    write_classification(result_tsv_lines, tsv_res_path)
    # Write accessions and names of hits to TSV file 'hits_to_download.tsv'
    write_hits_to_download(acc_dict, acc_fpath)

    # Update summary information
    seqs_processed[0] += len( packet["qual"] )
    pack_to_send[0] += 1
    remove_tmp_files(tmp_fpath)

# end def _handle_result


def retrieve_ready_job(saved_RID, packet, packet_size, packet_mode, pack_to_send, seqs_processed,
        fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
        blast_algorithm, user_email, organisms,
        acc_dict, out_of_n):
    """
    :param saved_RID: saved Request ID from previous run;
    :type saved_RID: str;
    :param packet: "packet" dictionary described in "barapost-prober.py" before the kernel loop:
    :type packet: dict;
    :param packet_size: size of the packet (see option '-c' for definition);
    :type packet_size: int;
    :param packet_mode: packet forming mode (see option '-c' for definition);
    :type packet_mode: int;
    :param pack_to_send: ordinal number of packet to send (it is list rather that in because it should be mutable);
    :type pack_to_send: list<int>;
    :param seqs_processed: nuber of sequnces processed (it is list rather that in because it should be mutable);
    :type seqs_processed: list<int>;
    :param fq_fa_path: path to current input file;
    :type fq_fa_path: str;
    :param tmp_fpath: path to current temporary file;
    :type tmp_fpath: str;
    :param taxonomy_path: path to taxonomt file;
    :type taxonomy_path: str;
    :param tsv_res_path: path to current classification file;
    :type tsv_res_path: str;
    :param acc_fpath: path to file 'hits_to_download.tsv';
    :type acc_fpath: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    :param blast_algorithm: BLAST algorithm to use (see option '-a');
    :type blast_algorithm: str;
    :param user_email: user email ot send with request;
    :type user_email: str;
    :param organisms: list of strings performing 'nt' database slices;
    :type organisms: list<str>;
    :param acc_dict: accession dictionary for writing to 'hits_to_download.tsv';
    :type acc_dict: dict<str: (str, int)>;
    :param out_of_n: dictionary for printing how many packets left;
    :type out_of_n: dict<str: str, str: int>;
    """


    resume_rtoe = 0 # we will not sleep at the very beginning of resumption

    # Get BLAST XML response
    align_xml_text, error = wait_for_align(saved_RID, resume_rtoe,
        pack_to_send, os.path.basename(fq_fa_path), logfile_path)

    if error.code == 0:
        # OK -- results are retrieved.
        # Get result tsv lines

        _handle_result(align_xml_text, packet, taxonomy_path,
            tsv_res_path, acc_dict, acc_fpath, seqs_processed, pack_to_send,
            tmp_fpath, logfile_path)

        return False
    elif error.code == 2:
        # If NCBI BLAST server rejects the request due to too large amount of data in it --
        #    split packet into two or, if there is only one sequence in it -- prune this sequence.
        # Then resend the request.

        _split_and_resubmit(packet, packet_size, packet_mode, pack_to_send, seqs_processed,
            fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
            blast_algorithm, user_email, organisms,
            acc_dict, out_of_n)
        return False
    else:
        return True
    # end if
# end def retrieve_ready_job


def submit(packet, packet_size, packet_mode, pack_to_send, seqs_processed,
        fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
        blast_algorithm, user_email, organisms,
        acc_dict, out_of_n):
    """
    :param packet: "packet" dictionary described in "barapost-prober.py" before the kernel loop:
    :type packet: dict;
    :param packet_size: size of the packet (see option '-c' for definition);
    :type packet_size: int;
    :param packet_mode: packet forming mode (see option '-c' for definition);
    :type packet_mode: int;
    :param pack_to_send: ordinal number of packet to send (it is list rather that in because it should be mutable);
    :type pack_to_send: list<int>;
    :param seqs_processed: nuber of sequnces processed (it is list rather that in because it should be mutable);
    :type seqs_processed: list<int>;
    :param fq_fa_path: path to current input file;
    :type fq_fa_path: str;
    :param tmp_fpath: path to current temporary file;
    :type tmp_fpath: str;
    :param taxonomy_path: path to taxonomt file;
    :type taxonomy_path: str;
    :param tsv_res_path: path to current classification file;
    :type tsv_res_path: str;
    :param acc_fpath: path to file 'hits_to_download.tsv';
    :type acc_fpath: str;
    :param logfile_path: path to log file;
    :type logfile_path: str;
    :param blast_algorithm: BLAST algorithm to use (see option '-a');
    :type blast_algorithm: str;
    :param user_email: user email ot send with request;
    :type user_email: str;
    :param organisms: list of strings performing 'nt' database slices;
    :type organisms: list<str>;
    :param acc_dict: accession dictionary for writing to 'hits_to_download.tsv';
    :type acc_dict: dict<str: (str, int)>;
    :param out_of_n: dictionary for printing how many packets left;
    :type out_of_n: dict<str: str, str: int>;
    """

    s_letter = 's' if len(packet["qual"]) != 1 else ''
    printl(logfile_path, "\nGoing to BLAST (" + blast_algorithm + ")")

    # Count base pairs in packet
    lines = filter(lambda x: not x.startswith('>'), packet["fasta"].splitlines())
    totalbp = len(''.join(map(lambda x: x.strip(), lines)))
    totalbp = "{:,}".format(totalbp)
    del lines

    printl(logfile_path, "Request number {}{}. Sending {} sequence{} ({} b.p. totally).".format(pack_to_send[0],
        out_of_n["msg"], len(packet["qual"]), s_letter, totalbp))

    error = BlastError(-1)

    while error.code != 0: # until successfull attempt

        request = configure_request(packet["fasta"], blast_algorithm, organisms, user_email) # get the request

        # Send the request and get BLAST XML response.
        # 'align_xml_text' will be None if an error occurs.
        align_xml_text, error = send_request(request, pack_to_send, packet_size, packet_mode,
            os.path.basename(fq_fa_path), tmp_fpath, logfile_path)

        if error.code == 0:
            # Write results and leave the loop
            _handle_result(align_xml_text, packet, taxonomy_path,
                tsv_res_path, acc_dict, acc_fpath, seqs_processed, pack_to_send,
                tmp_fpath, logfile_path)

        elif error.code == 2:
            # If NCBI BLAST server rejects the request due to too large amount of data in it --
            #    split packet into two or, if there is only one sequence in it -- prune this sequence.
            # Then resend the request.

            _split_and_resubmit(packet, packet_size, packet_mode, pack_to_send, seqs_processed,
                fq_fa_path, tmp_fpath, taxonomy_path, tsv_res_path, acc_fpath, logfile_path,
                blast_algorithm, user_email, organisms,
                acc_dict, out_of_n)

            error = BlastError(0) # _split_and_resubmit will process packet successfully
        # end if
    # end while
# end def submit