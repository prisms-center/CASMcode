import sklearn.linear_model

class SelectionWrapper(object):
    """ A wrapper function to handle selectors like RandomizedLasso,
    which cannot parse an "estimator" input but otherwise have
    the appropriate properties

    *** ONLY USE THIS IF YOU KNOW WHAT YOU ARE DOING ***
    """
    def __new__(cls, estimator, name=None, **kwargs):
        if hasattr(sklearn.linear_model, name):
          return getattr(sklearn.linear_model, name)(**kwargs)
        raise AttributeError("Could not find: " + name)

    def __init__(self):
        pass
