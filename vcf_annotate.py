  def split_info_field(self, line):
        """Method for splitting a VCF info field into a dict of its keys and values.
        Arguments:
            line (str): Single info entry to be split into a dictionary.
        Returns:
            info (dict): Dictionary containing each info field as a key 
                         with its corresponding value.
        """
        #Get key-value pairs. 
        pairs = line.split(';')
        info = {}
        for pair in pairs:
            field = pair.split('=')
            # Add key-value to dictionary 
            info[field[0]] = field[1]

        return info

