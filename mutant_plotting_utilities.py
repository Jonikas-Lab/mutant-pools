#! /usr/bin/env python
"""
Plotting utilities specifically for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import unittest
import glob
import math
# other packages
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties
# my modules
import mutant_analysis_classes
from general_utilities import unpickle
import plotting_utilities


################################################# Help functions ##################################################
# TODO should these data-extraction functions be in plotting_utilities or somewhere else?

def _extract_readcounts_from_dataset(dataset, perfect_only=False, total_and_perfect=False):
    """ Extract a readcount list from a mutant dataset. 

    Normally the returned list will contain total readcounts for each mutant;
     if perfect_only, the perfect instead of total readcounts are used; 
     if total_and_perfect, readcount_list will contain (total, perfect) tuples insteaad of single values.
    """
    if total_and_perfect:   return [(m.total_read_count, m.perfect_read_count) for m in dataset]
    elif perfect_only:      return [m.perfect_read_count for m in dataset]
    else:                   return [m.total_read_count for m in dataset]


def _get_readcount_data_if_needed(dataset, perfect_only=False, total_and_perfect=False):
    """ Return total/perfect readcount list from dataset object; if dataset is a number list already, just return it.

    For extraction details see _extract_readcounts_from_dataset.
    Of course, if the dataset is just a list of numbers, there's no way of checking whether those numbers are
     total or perfect-only, so the caller has to deal with that. 
    """
    # if dataset is a dataset object, extract the desired readcount list from it
    if isinstance(dataset, mutant_analysis_classes.Insertional_mutant_pool_dataset):
        return _extract_readcounts_from_dataset(dataset, perfect_only, total_and_perfect)
    # otherwise, check if dataset is a list of numbers (or number pairs if total_and_perfect): 
    #  if yes, return it, otherwise give an error.
    else:
        try:
            if total_and_perfect:
                assert all([len(x)==2 for x in dataset])
                sum(dataset[1]+dataset[2]+dataset[-1])
            else:
                sum(dataset[:10]+dataset[-10:])
            return dataset
        except TypeError, AssertionError:
            raise Exception("dataset format not recognized, in _get_readcount_data_if_needed!")


def get_all_datasets_glob(glob_pattern, split_filenames_on=None, readcounts_only=False, perfect_only=False, total_and_perfect=False):
    """ Get a name:dataset or name:readcount_list dict for all the files matching glob_pattern. 

    Files read using mutant_analysis_classes.read_mutant_file.

    By default returen a name:dataset dict, with full dataset objects. 
    If readcounts_only, return a name:readcount_list dict instead, to use less memory - 
        For extraction details see _extract_readcounts_from_dataset.

    By default dataset_name is file basename with no extension; if split_filenames_on is not None, 
     split the file basename on the value and use the first component as dataset_name.
    """
    all_datasets = {}
    for infile in glob.glob(glob_pattern):
        filename = splitext(os.path.basename(infile))[0]
        if split_filenames_on is not None:
            filename = filename.split(split_filenames_on)[0]
        dataset = mutant_analysis_classes.read_mutant_file(infile)
        if readcounts_only: all_datasets[filename] = _extract_readcounts_from_dataset(dataset, perfect_only, total_and_perfect)
        else:               all_datasets[filename] = dataset
    return all_datasets


############################################## Single-plot functions ##############################################


def plot_readcounts_sorted(dataset_name_list, color_dict=None, perfect_only=False, x_max=None, y_max=None, y_min=0.7, 
                           log_y=True, legend=True, legend_kwargs=None):
    """ Basic sorted readcount plots - multiple in a single plot, takes a (dataset,name) list. """
    for dataset,name in dataset_name_list:
        readcounts = _get_readcount_data_if_needed(dataset, perfect_only=perfect_only)
        kwargs = {} if color_dict is None else {'c': color_dict[name]}
        _ = mplt.plot(sorted(readcounts), '.', linewidth=0, label=name, markersize=2, **kwargs)
    mplt.xlabel('mutants sorted by readcount\n(independently for each dataset)')
    mplt.ylabel('readcount (%s)'%('log' if log_y else 'linear'))
    if log_y:   mplt.yscale('log')
    else:       mplt.yscale('linear')
    plotting_utilities.set_axes_limits(None, x_max, y_min, y_max)
    if legend:
        if legend_kwargs is None:   mplt.legend(loc=4, ncol=3, numpoints=3, handlelength=1.2, handletextpad=.7, columnspacing=1, 
                                                prop=FontProperties(size='medium'))
        else:                       mplt.legend(**legend_kwargs)
      

def plot_readcounts_cumulative(dataset_name_list, color_dict=None, linewidth_dict=None, linestyle_dict=None, perfect_only=False, 
                               x_min=None, x_max=None, y_max=None, y_min=None, log_x=True, legend=True, legend_kwargs=None):
    """ ROC-like "cumulative histogram" - multiple in a single plot, takes a (dataset,name) list. """
    for dataset,name in dataset_name_list:
        N_mutants = len(dataset)
        # adding a "0th mutant at readcount 1" data point to make sure all lines start at readcount 1, 
        #  even if the first real mutant has 100 reads - looks better that way, and more consistent with the nature of the plot.
        readcounts = [1] + sorted(_get_readcount_data_if_needed(dataset, perfect_only=perfect_only))
        mutant_percentages = [0] + [n/N_mutants*100 for n in range(N_mutants)]
        kwargs = {} 
        if color_dict is not None:      kwargs.update({'c': color_dict[name]})
        if linewidth_dict is not None:  kwargs.update({'linewidth': linewidth_dict[name]})
        if linestyle_dict is not None:  kwargs.update({'linestyle': linestyle_dict[name]})
        _ = mplt.plot(readcounts, mutant_percentages, label=name, **kwargs)
    mplt.xlabel('readcount (%s)'%('log' if log_x else 'linear'))
    if log_x:   mplt.xscale('log')
    else:       mplt.xscale('linear')
    mplt.ylabel('% of mutants with that readcount or less')
    plotting_utilities.set_axes_limits(x_min, x_max, y_min, y_max)
    if legend:
        if legend_kwargs is None:   mplt.legend(loc=4, ncol=3, numpoints=3, handlelength=1.2, handletextpad=.7, columnspacing=1, 
                                                prop=FontProperties(size='medium'))
        else:                       mplt.legend(**legend_kwargs)
  

def plot_readcounts_hist(dataset_name_list, color_dict=None, Nbins=100, histtype='bar', perfect_only=False, log_x=True, log_y=False, 
                         readcount_max=None, y_max=None, y_min=None, no_edges=False, legend=True, legend_kwargs=None):
    """ Normal histogram (any type), linear or logscale - multiple in a single plot, takes a (dataset,name) list. """
    for dataset,name in dataset_name_list:
        readcounts = _get_readcount_data_if_needed(dataset, perfect_only=perfect_only)
        if readcount_max:   readcounts = [c for c in readcounts if c<=readcount_max]
        if log_x:           readcounts = [math.log10(c) for c in readcounts]
        if log_x:           bin_range = (0, max(readcounts))
        else:               bin_range = (1, max(readcounts))
        kwargs = {} 
        if color_dict is not None:  kwargs.update({'facecolor': color_dict[name]})
        if no_edges:                kwargs.update({'edgecolor': 'none'})
        # using range to make sure that if we have one dataset with some 1-read mutants and one with only 10+-read mutants, 
        #  they both get the same bins, rather than N bins stretched over the (1,max) and (10,max) ranges respectively.
        _ = mplt.hist(readcounts, histtype=histtype, linewidth=0.3, bins=Nbins, range=bin_range, label=name, log=log_y, 
                      align='mid', **kwargs)
    mplt.xlabel('readcount (%s) (range binned into %s bins)'%(('log' if log_x else 'linear'), Nbins))
    mplt.ylabel('number of mutants with that readcount (%s)'%('log' if log_y else 'linear'))
    # TODO if x_log, change the xticklabels (and maybe add ticks?) to look like log!
    plotting_utilities.set_axes_limits(None, readcount_max, y_min, y_max)
    if legend:
        mplt.legend(loc=1, ncol=3, numpoints=3, handlelength=1.2, handletextpad=.7, columnspacing=1, 
                    prop=FontProperties(size='medium'))
  

############################################# Multi-plot functions ##############################################


# this is OLD, probably doesn't work - MAYBE-TODO make it work if I actually want to?
def plot_sample_row_multi_plots(sample_name, sample_N, total_samples, first_cumulative=True, 
                                xmax_ymax_Nbins_logx_logy_list=[(None,None,100,True,False),(None,None,100,False,False)], 
                                histtype='bar', if_xlabels=True):
    """ Plot one sample readcount info with multiple methods: cumulative histogram, and any number of normal histograms with different settings. """
    N_plots = len(xmax_ymax_Nbins_logx_logy_list) + int(first_cumulative)
    curr_plot = sample_N * N_plots + 1
    if first_cumulative:
        mplt.subplot(total_samples, N_plots, curr_plot)
        plot_data_cumulative([sample_name], legend=False)
        if not if_xlabels:
            mplt.xlabel('')
            plotting_utilities.remove_xticklabels()
        mplt.ylabel(sample_name)
        curr_plot += 1
    for xmax,ymax,Nbins,log_x,log_y in xmax_ymax_Nbins_logx_logy_list:
        mplt.subplot(total_samples, N_plots, curr_plot)
        plot_data_hist([sample_name], Nbins=Nbins, logx=log_x, logy=log_y, maxcount=xmax, histtype=histtype, legend=False)
        if ymax is not None:
            mplt.ylim((0,ymax))
        if not if_xlabels:
            mplt.xlabel('')
            plotting_utilities.remove_xticklabels()
        if curr_plot != sample_N * N_plots + 1:
            mplt.ylabel('')
            plotting_utilities.remove_yticklabels()
        else:
            mplt.ylabel(sample_name)
        curr_plot += 1


################################################# Testing etc ##################################################

class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__(self):
        sys.exit("NO UNIT-TESTS FOR THIS MODULE")
    # LATER-TODO add unit-tests!


if __name__=='__main__':
    """ If module is run directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
