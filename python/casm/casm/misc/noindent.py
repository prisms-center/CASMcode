import json
import six
import uuid

# ---------------------------------------------------
# Some code to keep parts of a json string from being spread across multiple lines
# code from: http://stackoverflow.com/questions/13249415/can-i-implement-custom-indentation-for-pretty-printing-in-python-s-json-module
# answer: http://stackoverflow.com/a/25935321
# code by: http://stackoverflow.com/users/247623/erik-allik

class NoIndent(object):
    def __init__(self, value):
        self.value = value


class NoIndentEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        super(NoIndentEncoder, self).__init__(*args, **kwargs)
        self.kwargs = dict(kwargs)
        del self.kwargs['indent']
        self._replacement_map = {}

    def default(self, o):
        if isinstance(o, NoIndent):
            key = uuid.uuid4().hex
            self._replacement_map[key] = json.dumps(o.value, **self.kwargs)
            return "@@%s@@" % (key,)
        else:
            return super(NoIndentEncoder, self).default(o)

    def encode(self, o):
        result = super(NoIndentEncoder, self).encode(o)
        for k, v in six.iteritems(self._replacement_map):
            result = result.replace('"@@%s@@"' % (k,), v)
        return result
# ---------------------------------------------------
