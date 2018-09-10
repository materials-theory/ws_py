
class OutcarHandler(object):
    # TODO : substitute all functions using regular expression

    @staticmethod
    def param_from_outcar(keyword):
        with open('OUTCAR', 'r') as out:
            for line in out:
                if keyword in line:
                    return line.split()[2]
