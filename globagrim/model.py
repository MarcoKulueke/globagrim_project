from .constants import global_const
from . import globagrim
from . import visualization


def run(**kwargs):

    global_const.IEXP = kwargs.get("IEXP", global_const.IEXP)

    global_const.NJ = kwargs.get("NJ", global_const.NJ)
    global_const.NK = kwargs.get("NK", global_const.NK)
    global_const.NL = kwargs.get("NL", global_const.NL)

    global_const.DT = kwargs.get("DT", global_const.DT)
    global_const.TF = kwargs.get("TF", global_const.TF)

    globagrim.globagrim()


def plot(**kwargs):
    
    var_name = kwargs.get("var_name", 'PSG')
    time_step = kwargs.get('time_step', 0)
    
    visualization.plot(var_name, time_step)
