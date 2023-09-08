#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample_col="sample",
        ref_col="reference",
        index_col="index",
        first_col=False,
        second_col=False,
        bam_col = False,
        single_col=False,
        whitelist_col=False,
        cellnumber_col=False,
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._bam_col = bam_col
        self._single_col = single_col
        self._ref_col = ref_col
        self._index_col = index_col
        self._whitelist_col = whitelist_col
        self._cellnumber_col = cellnumber_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        if self._first_col != False:
            self._validate_first(row)
        if self._second_col != False:
            self._validate_second(row)
            self._validate_pair(row)
        if self._whitelist_col != False:
            self._validate_whitelist_cellnumber(row)
        if self._bam_col != False:
            self._validate_bam_format(row)
        # TODO check that reference is fasta or fa file and either reference or index is given if mode is not bulk
        # TODO check that combination of reference and index are consistent
        if self._ref_col != False:
            pass

        self._seen.add((row[self._sample_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("At least the first FASTQ file is required.")
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        if len(row[self._second_col]) > 0:
            self._validate_fastq_format(row[self._second_col])

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            # row[self._single_col] = False
            first_col_suffix = Path(row[self._first_col]).suffixes[-2:]
            second_col_suffix = Path(row[self._second_col]).suffixes[-2:]
            if first_col_suffix != second_col_suffix:
                raise AssertionError("FASTQ pairs must have the same file extensions.")
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FORMATS):
            raise AssertionError(
                f"The FASTQ file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )
    def _validate_bam_format(self, row):
        filename = row[self._bam_col]
        if not filename.endswith(".bam"):
            raise AssertionError(
                f"The BAM file has an unrecognized extension: {filename}\n"
                f"It should be .bam"
            )
        
    def _validate_whitelist_cellnumber(self, row):
        """Assert that either a whitelist or the expected number of cells are given."""
        if (len(row[self._whitelist_col]) <= 0 and len(row[self._cellnumber_col]) <= 0) or (len(row[self._whitelist_col]) > 0 and len(row[self._cellnumber_col]) > 0):
            raise AssertionError("Either whitelist or cellnumber must be provided")
        # TODO check that cellnumber is integer
        # TODO check that whitelist has ending tsv txt or csv

    def _validate_ref_index(self, row):
        pass

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename all samples to have a suffix of _T{n}, where n is the
        number of times the same sample exist, but with different FASTQ files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The sample name must be unique.")
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            seen[sample] += 1


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical("The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(mode, file_in, file_out, input_type, pipeline, reference):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            sample,fastq_1,fastq_2
            SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz
            SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz
            SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,

    .. _viral recon samplesheet:
        https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

    """
    if mode == "single-bulk":
        if reference:
            required_columns = {"sample", "fastq_1", "reference", "index"}
        else:
            required_columns = {"sample", "fastq_1"}

    elif mode == "paired-bulk": 
        if reference:
            required_columns = {"sample", "fastq_1", "fastq_2", "reference", "index"}
        else:
            required_columns = {"sample", "fastq_1", "fastq_2"}

    elif mode == "single-cell": 
        if input_type == "fastq":
            if pipeline == "saw":
                required_columns = {"sample", "fastq_1", "reference", "index"}
            else:
                required_columns = {"sample", "fastq_1", "fastq_2", "reference", "index", "whitelist", "cellnumber"}

        elif input_type == "bam":
            required_columns = {"sample", "bam", "reference", "index"}
        else:
            logger.critical(f"Unrecognized parameter combination: {mode}, {input_type}, {pipeline}.")
            sys.exit(1)
    else:
        logger.critical(f"Unrecognized parameter combination: {mode}, {input_type}, {pipeline}.")
        sys.exit(1)

    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        # Initialize the RowChecker with the columns depending on which mode the pipeline is run in.
        if mode == "single-bulk":
            if reference:
                checker = RowChecker(first_col="fastq_1")
            else:
                RowChecker(first_col="fastq_1", ref_col=False, index_col=False)

        elif mode == "paired-bulk":
            if reference:
                checker = RowChecker(first_col="fastq_1", second_col="fastq_2")
            else:
                RowChecker(first_col="fastq_1", ref_col=False, index_col=False)
            
        elif mode == "single-cell": 
            if input_type == "fastq":
                if pipeline == "saw":
                    checker = RowChecker(first_col="fastq_1")
                else:
                    checker = RowChecker(first_col="fastq_1", second_col="fastq_2", whitelist_col="whitelist", cellnumber_col="cellnumber")

            elif input_type == "bam":
                checker = RowChecker(bam_col="bam")

        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
    checker.validate_unique_samples()
    header = list(reader.fieldnames)
    # header.insert(1, "single_end")
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "mode",
        metavar="MODE",
        type=str,
        help="Mode pipeline is run in.",
    )
    parser.add_argument(
        "reference",
        metavar="REFERENCE",
        type=bool,
        help="Whether to align barcodes to a reference.",
    )
    parser.add_argument(
        "-t",
        "--input_type",
        metavar="INPUT_TYPE",
        required=False,
        type=str,
        help="Whether input is in fastq or bam format.",
        choices=["fastq", "bam"]
    )
    parser.add_argument(
        "-p",
        "--pipeline",
        metavar="PIPELINE",
        required=False,
        type=str,
        help="Pipeline with which input was generated (saw, cellranger or starsolo).",
        choices=["saw", "cellranger", "starsolo"]
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(mode=args.mode, file_in=args.file_in, file_out=args.file_out, input_type=args.input_type, pipeline=args.pipeline, reference=args.reference)


if __name__ == "__main__":
    sys.exit(main())
