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
import general_utilities
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


############################################## Plotting readcounts ##############################################

######### Single-plot functions

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
    plotting_utilities.remove_half_frame()
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
    plotting_utilities.remove_half_frame()
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
    plotting_utilities.remove_half_frame()
    if legend:
        mplt.legend(loc=1, ncol=3, numpoints=3, handlelength=1.2, handletextpad=.7, columnspacing=1, 
                    prop=FontProperties(size='medium'))
  

######### Multi-plot functions

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

######################################## Plotting adjacent mutant data ###########################################

def _get_singles_from_counter(val_count_dict, max_val):
    """ Given a val:count dict, return a val list with each val present count times. (for creating histograms). """
    val_singles_list = []
    for val,count in val_count_dict.items():
        if max_val is None or val <= max_val:
            val_singles_list.extend([val]*count)
    return val_singles_list


colors_by_adjacent_category = {'adjacent-same-strand': 'red', 'adjacent-opposite-both': 'cyan', 'same-pos-opposite': 'orange', 
                               'adjacent-opposite-toward': 'green', 'adjacent-opposite-away': 'blue' }


def plot_adjacent_distance_histogram(dataset, incl_same_strand=True, incl_opposite_both=True, incl_opposite_separate=True, 
                                     incl_same_pos_opposite=False, max_distance=None, N_bins=100, symlog_y=False, symthresh=100):
    """ Step-histogram of the number of adjacent mutant pairs by distance, separated by type. 

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, which had its count_adjacent_mutants 
     method ran at some point to fill the relevant summary fields with the information for this plot.
    All the incl_* variables are True/False and govern which adjacent mutant types should be included. 
     (Be careful with incl_same_pos_opposite - those are always distance 0, so it's hard to show them meaningfully in comparison
       to the ones with a wider ditance range unless each distance has its own bin, like with max_distance 100 and N_bins 101). 
    Only include mutant pairs with distance at most max_distance, if given.
    If symlog_y, the y axis will be symlog-scale (linear up to symthresh, exponential afterward). 
    """
    datasets, labels, colors = [], [], []
    if incl_same_strand:    
        datasets.append(_get_singles_from_counter(dataset.summary.adjacent_same_strand_dict, max_distance))
        labels.append('adjacent-same-strand')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_opposite_separate:
        datasets.append(_get_singles_from_counter(dataset.summary.adjacent_opposite_toward_dict, max_distance))
        labels.append('adjacent-opposite-toward')
        colors.append(colors_by_adjacent_category[labels[-1]])
        datasets.append(_get_singles_from_counter(dataset.summary.adjacent_opposite_away_dict, max_distance))
        labels.append('adjacent-opposite-away')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_opposite_both:
        adjacent_opposite_both = general_utilities.add_dicts_of_ints(dataset.summary.adjacent_opposite_toward_dict, 
                                                                     dataset.summary.adjacent_opposite_away_dict)
        datasets.append(_get_singles_from_counter(adjacent_opposite_both, max_distance))
        labels.append('adjacent-opposite-both')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_same_pos_opposite:
        datasets.append([0]*dataset.summary.same_position_opposite)
        labels.append('same-pos-opposite')
        colors.append(colors_by_adjacent_category[labels[-1]])
        if max_distance>=N_bins:
            print("Warning: if each distance doesn't get its own bin (N_bins = max_distance+1), "
                  +"the same-pos-opposite count won't be fairly represented!")
    xmin = min(sum(datasets, []))
    xmax = max(sum(datasets, []))
    xspread = xmax - xmin
    # have to do bins by hand - TODO or do I?
    #bins = [(low + i/N_bins*spread) for i in range(N_bins+1)]
    hist_data, bin_edges, hist_patches = mplt.hist(datasets, label=labels, color=colors, bins=N_bins, histtype='step')
    plotting_utilities.remove_half_frame()
    mplt.xlabel('distance between the two mutants in a pair (bp)' + (' (binned into %s bins)'%N_bins if (xspread+1)>N_bins else ''))
    mplt.ylabel('number of mutant pairs with given distance')
    # sometimes mplt.hist gets the axis limits wrong, so fix them by hand
    ymax = max(sum([list(d) for d in hist_data], []))
    mplt.xlim(xmin-xspread/N_bins, xmax)
    mplt.ylim(0-ymax/100, ymax*1.05)
    mplt.legend(title='pair categories by relative orientation:', prop=FontProperties(size='medium'))


def plot_adjacent_dist_1_bars(dataset, incl_same_strand=True, incl_opposite_both=True, incl_opposite_separate=True, 
                              incl_same_pos_opposite=True, logscale_x=False):
    """ Horizontal bar-plot of the number of 0-1bp adjacent mutant pairs, separated by type. 

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, which had its count_adjacent_mutants 
     method ran at some point to fill the relevant summary fields with the information for this plot.
    All the incl_* variables are True/False and govern which adjacent mutant types should be included. 
    If logscale_x, the x axis will be log-scale. 
    """
    counts, labels, colors = [], [], []
    if incl_same_strand:    
        counts.append(dataset.summary.adjacent_same_strand_dict[1])
        labels.append('adjacent-same-strand')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_opposite_both:
        counts.append(dataset.summary.adjacent_opposite_toward_dict[1] + dataset.summary.adjacent_opposite_away_dict[1])
        labels.append('adjacent-opposite-both')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_opposite_separate:
        counts.append(dataset.summary.adjacent_opposite_toward_dict[1])
        labels.append('adjacent-opposite-toward')
        colors.append(colors_by_adjacent_category[labels[-1]])
        counts.append(dataset.summary.adjacent_opposite_away_dict[1])
        labels.append('adjacent-opposite-away')
        colors.append(colors_by_adjacent_category[labels[-1]])
    if incl_same_pos_opposite:
        counts.append(dataset.summary.same_position_opposite)
        labels.append('same-pos-opposite')
        colors.append(colors_by_adjacent_category[labels[-1]])
    # use everything reversed, so that the first one is on top instead of bottom (convert to list, matplotlib chokes on iterators)
    mplt.barh(range(len(counts)), list(reversed(counts)), color=list(reversed(colors)), align='center', log=logscale_x)
    mplt.yticks(range(len(labels)), list(reversed(labels)))
    mplt.xlabel('number of mutant pairs in given category')
    plotting_utilities.remove_half_frame()


def _get_sorted_ratios_from_dict(ratio_distance_dict, min_distance, max_distance):
    """ Given a distance:readcount_pair_list dict, return a list of readcount_pair ratios for distances between min and max. 
    """
    ratio_list = []
    for dist in ratio_distance_dict.keys():
        if min_distance <= dist <= max_distance:
            for x,y in ratio_distance_dict[dist]:
                x,y = sorted([x,y])
                ratio_list.append(y/x)
    return sorted(ratio_list)


def plot_adjacent_readcount_ratios(dataset, distance_cutoffs, distance_linestyles=None, incl_same_strand=True, 
                           incl_opposite_both=True, incl_opposite_separate=True, incl_same_pos_opposite=True, logscale_y=True):
    """ Plot sorted readcount ratios of adjacent mutant pairs, separated by type and optionally distance categories.

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, which had its count_adjacent_mutants 
     method ran at some point to fill the relevant summary fields with the information for this plot.
    Distance_cutoffs should be a list of numbers: [1,10,1000] will result in plotting separate lines for mutant pairs with 
     distance 1, 2-10, and 11-1000 for each type (except same-pos-opposite, since those always have a distance of 0)
    Distance_linestyles should be a list of valid matplotlib.pyplot.plot linestyle values, one for each cutoff in distance_cutoffs.
    All the incl_* variables are True/False and govern which adjacent mutant types should be included. 
    If logscale_y, the y axis will be log-scale. 
    """
    # LATER-TODO do I actually want sorted readcount ratios, or would a histogram or something be better?  But histograms don't have linestyles... I could just make histogram-like data and use mplt.plot to plot it, though.
    if distance_linestyles is None:
        if len(distance_cutoffs)==2:    distance_linestyles = ['-', ':']
        elif len(distance_cutoffs)==3:  distance_linestyles = ['-', '--', ':']
        elif len(distance_cutoffs)==4:  distance_linestyles = ['-', '--', '-.', ':']
        else:   raise Exception("If there aren't 2-4 distance_cutoffs, need to provide distance_linestyles!")
    ratio_lists, labels, colors, linestyles = [], [], [], []
    distance_mins = [0] + [x+1 for x in distance_cutoffs[:-1]]
    if incl_same_strand:    
        category = 'adjacent-same-strand'
        for distance_min,distance_max,linestyle in zip(distance_mins, distance_cutoffs, distance_linestyles):
            ratio_lists.append(_get_sorted_ratios_from_dict(dataset.summary.adjacent_same_strand_readcounts_dict, 
                                                            distance_min, distance_max))
            labels.append('%s, dist %s-%s'%(category, distance_min, distance_max))
            colors.append(colors_by_adjacent_category[category])
            linestyles.append(linestyle)
    if incl_opposite_both:
        category = 'adjacent-opposite-both'
        for distance_min,distance_max,linestyle in zip(distance_mins, distance_cutoffs, distance_linestyles):
            ratio_lists.append(sorted(sum([_get_sorted_ratios_from_dict(D, distance_min, distance_max) 
                                           for D in (dataset.summary.adjacent_opposite_toward_readcounts_dict, 
                                                     dataset.summary.adjacent_opposite_away_readcounts_dict,)], [])))
            labels.append('%s, dist %s-%s'%(category, distance_min, distance_max))
            colors.append(colors_by_adjacent_category[category])
            linestyles.append(linestyle)
    if incl_opposite_separate:
        category = 'adjacent-opposite-toward'
        for distance_min,distance_max,linestyle in zip(distance_mins, distance_cutoffs, distance_linestyles):
            ratio_lists.append(_get_sorted_ratios_from_dict(dataset.summary.adjacent_opposite_toward_readcounts_dict, 
                                                            distance_min, distance_max))
            labels.append('%s, dist %s-%s'%(category, distance_min, distance_max))
            colors.append(colors_by_adjacent_category[category])
            linestyles.append(linestyle)
        category = 'adjacent-opposite-away'
        for distance_min,distance_max,linestyle in zip(distance_mins, distance_cutoffs, distance_linestyles):
            ratio_lists.append(_get_sorted_ratios_from_dict(dataset.summary.adjacent_opposite_away_readcounts_dict, 
                                                            distance_min, distance_max))
            labels.append('%s, dist %s-%s'%(category, distance_min, distance_max))
            colors.append(colors_by_adjacent_category[category])
            linestyles.append(linestyle)
    if incl_same_pos_opposite:
        category = 'same-pos-opposite'
        ratio_lists.append(_get_sorted_ratios_from_dict({0: dataset.summary.same_position_opposite_readcounts}, 0, 0))
        labels.append('%s (always dist 0)'%(category))
        colors.append(colors_by_adjacent_category[category])
        linestyles.append(distance_linestyles[0])
    for ratio_list,label,color,linestyle in zip(ratio_lists, labels, colors, linestyles):
       mplt.plot([(y+1)/len(ratio_list) for y in range(len(ratio_list))], ratio_list, 
                 color=color, linestyle=linestyle, label='%s - %s pairs'%(label, len(ratio_list)))
    if logscale_y:   mplt.yscale('log')
    plotting_utilities.remove_half_frame()
    mplt.ylabel('readcount ratio between the two mutants in a pair')
    mplt.xlabel('all adjacent mutant pairs, sorted by readcount ratio, normalized to 100%')
    smallfont = FontProperties(size='smaller')
    mplt.legend(prop=smallfont, loc=2)
    mplt.title('Ratios between the readcounts of adjacent mutant pairs,\nby category and distance')


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
