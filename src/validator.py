class CustomError(BaseException):

    def __init__(self,msg=None):
        super().__init__(msg)
        

def validate_input_data(df):
    '''
    This method validates the input dataframe for required columns and missing values.

    :param df: This is the input dataframe to be validated
    '''

    # Validate for required columns - 'protein' column is required and at least one of 'motif' or 'site' column is required.
    column_names = df.columns.to_list()
    if 'protein' not in column_names:
        raise CustomError(msg="Input file error:: Column, 'protein' is missing")
    elif (('motif' not in column_names) and ('site' not in column_names)):
        raise CustomError(msg="Input file error:: Either column - 'motif' or 'site' must be present")
    
    # expected_col_names = ['protein','motif', 'site', 'logFC','p-value']
    # if len(column_names) > 2:
    #     for col_name in column_names:
    #         if col_name not in expected_col_names:
    #              raise CustomError(msg=f"Input file error:: Incorrect column {col_name}")

    # If protein column has missing values, raise error
    col = 'protein'
    if df[col].isnull().any():
        raise CustomError(msg=f"Input file error:: Missing values/None values found in column, '{col}'")
    
    # If site column is present, validate for missing values and correct format of site (S/T/Y followed by a number)
    if 'site' in df.columns:
        if df['site'].isnull().any():
            raise CustomError(msg=f"Input file error:: Missing values/None values found in column, 'site'")
        sites = df['site'].to_list()
        for site in sites:
            if site[0].upper() not in ['S','T','Y']:
                raise CustomError(msg=f"Input file error:: Site, '{site}' does not start with S/T/Y")
            if not site[1:].isdigit():
                raise CustomError(msg=f"Input file error:: Site, '{site}' does not have a valid number after the residue")
    else:
        raise CustomError(msg=f"Input file error:: Column, 'site' is missing.")
        # If site column is not present, then motif column must be present. Validate for missing values and correct format of motif (length of 9 with S/T/Y at the center)
        # if 'motif' in df.columns:
        #     if df['motif'].isnull().any():
        #         raise CustomError(msg=f"Input file error:: Missing values/None values found in column, 'motif'")
        #     motifs = df['motif'].to_list()
        #     for motif in motifs:
        #         if len(motif) != 9:
        #             raise CustomError(msg=f"Input file error:: Motif, '{motif}' is not of length 9")
        #         if motif[4] not in ['S','T','Y']:
        #             raise CustomError(msg=f"Input file error:: Motif, '{motif}' does not have S/T/Y at the center position")

    # If logFC column is present, validate for missing values
    if 'logFC' in df.columns:
        if df['logFC'].isnull().any():
            raise CustomError(msg=f"Input file error:: Missing values/None values found in column, 'logFC'")
    
    # If p-value column is present, validate for missing values
    if 'p-value' in df.columns:
        if df['p-value'].isnull().any():
            raise CustomError(msg=f"Input file error:: Missing values/None values found in column, 'p-value'")