class VCFMetadata:
    """Class for parsing metadata from a VCF file.
    Arguments:
        vcf_file (str): Path to VCF file passed on object creation.
    Attributes:
        file_name (str): Path to VCF file.
        raw_meta_data (list(str)): Raw lines containing metadata.
        headers (list(str)): List of header strings.
        meta_data_length (int): Number of lines containing metadata. 
    Raises:
        ValueError: Exception raised when no metadata can be parsed from
                    the input VCF file.
    """

    def __init__(self, vcf_file):
        self.file_name = vcf_file
        with open(vcf_file) as vcf_handler:
            # Find all meta data lines starting with #.
            self.raw_meta_data = [line.rstrip('\n') for line in vcf_handler
                                  if line.startswith('#')]
        # Get number of meta data lines.
        self.meta_data_length = len(self.raw_meta_data)
        # Raise ValueError if no meta data was able to be parsed.
        if self.meta_data_length == 0:
            raise ValueError('Failed to parse metadata from VCF file.\n'\
                             'Check input for proper formatting.')
        # Strip # and split last line of meta data to get headers.
        self.headers = self.raw_meta_data[-1][1:].split('\t')

        #TODO parse rest of meta data?
