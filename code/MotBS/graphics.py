# -*- coding: utf-8 -*-
"""

Created on Wed Aug 05 2011
Отделение панелей хвостового отсека: построение графиков 
@author: Юдинцев В. В.

"""

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pylab import *

from matplotlib import rc
rc('font',**{'family':'serif','size':'14'})
rc('text', usetex=True)
rc('text.latex',unicode=True)
rc('text.latex',preamble='\usepackage[utf8]{inputenc}')
rc('text.latex',preamble='\usepackage[russian]{babel}')
#rc('text.latex',preamble='\usepackage{amsmath}')
#rc('text.latex',preamble='\usepackage{amssymb}')

def figsize(wcm,hcm): 
    result = figure(figsize=(wcm/2.54,hcm/2.54))
    return result

def create_plot() :
    figure = figsize(18,12)    
    ax = figure.add_subplot(111)
    figure.subplots_adjust(top=0.95, bottom=0.05, left=0.08, right=0.94)    
    return figure, ax 
    
def label_gost(ax, xlabel, ylabel) :
    # Определяем интервал отображения для оси x
    vxi = ax.xaxis.get_view_interval()    
    # результатом является список vi[0] = xmin, vi[1] = xmax    
    
    def fmtx(x, pos=None):
        # если пришло время рисовать самое правое число на оси x
        # заменяем его на подпись к оси  
        if x == vxi[1] :
            return xlabel            
        # самое левое число оси не рисуем, т.к. оно слишком близко 
        # находится от подписи к оси ординат
        elif pos == 0 :
            return ''
        # в прочих случаях форматируем число на оси абсцисс
        else :
            return '%.3g' % x
        
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmtx))
    
    # Аналогичный код для оси y
    vyi = ax.yaxis.get_view_interval()    
    def fmty(y, pos=None):
        # print 'y:', x, pos, lymax
        if y == vyi[1] :            
            return ylabel            
        return '%.4g' % y
        
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmty))    
    
    # TODO: форматирование в одной функции
    


# =============================================================================
# 
# =============================================================================
   
path='/home/elim/Programs/MotBS/'        


stat   = np.genfromtxt(path+"result.csv", delimiter=",", autostrip=True)
stat_t = np.transpose(stat)

RTOD = 180.0/pi

fig, ax = create_plot()
ax.plot(stat_t[0],stat_t[1]*RTOD,'-',stat_t[0],stat_t[2]*RTOD,'--',stat_t[0],stat_t[3]*RTOD,'-.',linewidth=2)
label_gost(ax,'t, c','$\\varphi$, $^\circ$')
legend(('$\\varphi_1$', '$\\varphi_2$','$\\varphi_3$'))
grid(True)
savefig(path+'phi.png', dpi = 300)

fig, ax = create_plot()
ax.plot(stat_t[0],stat_t[4]*RTOD,'-',stat_t[0],stat_t[5]*RTOD,'--',stat_t[0],stat_t[6]*RTOD,'-.',linewidth=2)
label_gost(ax,'t, c','$\dot \\varphi$, $^\circ/c$')
legend(('$\dot \\varphi_1$', '$\dot \\varphi_2$','$\dot \\varphi_3$'))
grid(True)
savefig(path+'dphi.png', dpi = 300)

