def list_sample_models():
    """
    Returns all included sample models
    """
    return(['Bacillus-subtilis-168-WT',])
def get_sample_model(model):
    """
        Returns path to sample model

        Parameters
        ----------
        model: Name of model (folder with RBA-xml files in sample_models directory)
        
        Returns
        -------
        path to sample model (relative to jupyter notebook directory).
    """

    return('sample_models/{}'.format(model))