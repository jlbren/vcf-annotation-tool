from vcf_meta import VCFMetadata
import pandas

class VCFParser:
    """Class for parsing VCF files. Creates a pandas dataframe 
       and VCFMetadata object.
    Arguments:
        vcf_file (str): Path to VCF file passed on object creation.
    Attributes:
        meta_data: VCFMetadata object generated from input file.
        vcf_df: Pandas dataframe created from VCF data. 
    """

    def __init__(self, vcf_file):
        # Create VCFMetadata object. 
        self.meta_data = VCFMetadata(vcf_file)
        # Create a pandas dataframe from VCF data.
        # Use headers parsed by VCFMetadata, and skip metadata lines. 
        self.vcf_df = pandas.read_csv(self.meta_data.file_name,
                                      sep = '\t',
                                      names = self.meta_data.headers,
                                      skiprows = self.meta_data.meta_data_length)

  
