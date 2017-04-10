#!/usr/bin/env python

import numpy as np



def chi2(off_set,template_func,template_err_func, spec):
    """
    The module is to create chi^2.
    """
    shift = 1. + off_set/299792.458 # (1+v/c)
    chi_squared = np.sum((spec.flux[i] - template_func(wave*shift))**2/\
                             (spec.error[i]**2 + template_err_func(wave*shift)**2)\
                             for i,wave in enumerate(spec.wave)\
                             if spec.mask[i])
    return chi_squared
# If you want a reduced chi2 (chi2/dof), change to
#    appr_size = np.count_nonzero(spec.mask)
#    return chi_squared/(appr_size-1.0)


def detail_wl_array(old_array,resol):
    """
    The module is to create a new, with smaller step, wavelength array.
    
    (resol-1) - number of subpixels between the old pixels.
    """
    new_list = []
    len_old = len(old_array)
    resol_rev = 1./resol
    for i in range(len_old):
        new_list.append(old_array[i])
        if i != len_old-1:
            for j in range(1,resol):
                new_list.append(old_array[i] + resol_rev*j*(old_array[i+1] - old_array[i]))
    return np.array(new_list)


def off_set_plot(shift_list):
    """
    The module is to create a bar chart with distribution of the off-sets.
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib import rc
    #rc('font', family='Comic Sans MS')
    rc('axes', labelsize=16)
#    with open(shift_list) as f:# in case list is path to file voffset_result.dat
#        data = f.read().split('\n')
#    data = map(float,[data_i.split(' ')[1] for data_i in data])
    data = shift_list
    bin_min = np.min(data)
    bin_max = np.max(data)
    bin_step = 0.1
    if bin_min >= 0:
        bins = np.arange(bin_step*(bin_min//bin_step - 1), bin_step*(bin_max//bin_step + 2),bin_step)
    elif bin_max < 0:
        bins = np.arange(bin_step*(bin_min//bin_step), bin_step*(bin_max//bin_step + 3),bin_step)
    elif bin_min < 0 and bin_max >= 0:
        bins = np.arange(bin_step*(bin_min//bin_step), bin_step*(bin_max//bin_step + 2),bin_step)
    else:
        print "Something wrong with min/max bins!"
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    n, bins, patches = ax.hist(data, bins, histtype='bar', facecolor='#a6cee3',edgecolor = '#1f78b4',alpha=1)# rwidth=0.95
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.set_ylabel('Number of spectra')
    ax.set_xlabel('Velocity off-set, km/s')
    ax.set_xlim([np.min(bins)-0.95*bin_step,np.max(bins)+0.95*bin_step])
    ax.set_ylim([0.01,np.max(n)+0.3])
    fig.savefig('voffset_dist.pdf',bbox_inches='tight',pad_inches=0)
    print "................................................................................."
    print ""
    print "The file 'voffset_dist.pdf' has been created!"
    plt.close()
