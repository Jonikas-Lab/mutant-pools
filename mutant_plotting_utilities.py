#! /usr/bin/env python
"""
Plotting utilities specifically for mutant datasets and related things.  Module - running it directly just runs tests.
 -- Weronika Patena, 2012
"""

# standard library
from __future__ import division
import unittest
import glob
import os
from collections import defaultdict
import math
# other packages
import numpy
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties
# my modules
import mutant_analysis_classes
import general_utilities
import seq_basic_utilities
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

def readcounts_sorted_plot(dataset_name_list, color_dict=None, perfect_only=False, x_max=None, y_max=None, y_min=0.7, 
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
      

def readcounts_cumulative_plot(dataset_name_list, color_dict=None, linewidth_dict=None, linestyle_dict=None, perfect_only=False, 
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
  

def readcounts_histogram(dataset_name_list, color_dict=None, Nbins=100, histtype='bar', perfect_only=False, log_x=True, log_y=False, 
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
def sample_row_multi_plots(sample_name, sample_N, total_samples, first_cumulative=True, 
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


########################################## Plotting mutant positions #############################################

def _get_mutant_positions_from_dataset(dataset, strand='both'):
    """  Return chromosome_name:mutant_position_list for dataset.

    Dataset must be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance.
    Use the (known or assumed) position of the base before the insertion (min_position).
    If strand is 'both', take mutants regardless of strand; 
     if it's '+' or '-', take only mutants on that strand; all other values are illegal.
    """
    chromosome_position_dict = defaultdict(list)
    for mutant in dataset:
        position = mutant.position
        if strand=='both' or position.strand==strand:
            chromosome_position_dict[position.chromosome].append(position.min_position)
    return chromosome_position_dict


def _get_chromosome_lengths(genome_file=None):
    """ Return chromosome:length dictionary based on reading a genome fasta file. """
    if genome_file is None:
        genome_file = os.path.expanduser('~/experiments/reference_data/genomes_and_indexes/Chlre4-nm_chl-mit_cassette-pMJ013b.fa')
    chromosome_lengths = defaultdict(int)
    for header,seq in seq_basic_utilities.parse_fasta(genome_file):
      chromosome_lengths[header] = len(seq)
    return dict(chromosome_lengths)


def _get_chromosomes_by_type(all_chromosomes, include_scaffolds, include_cassette, include_other, output_dict=False):
    """ Filter all_chromosomes based on the include_* arguments; return a list of only ones of the desired types. 

    If output_dict is True, output a type:chromosome_list dictionary instead of a single chromosome_list.
    """
    chosen_chromosomes = defaultdict(list)
    for chromosome in all_chromosomes:
        chr_type = seq_basic_utilities.chromosome_type(chromosome)
        if (chr_type=='chromosome' or (chr_type=='scaffold' and include_scaffolds) or (chr_type=='cassette' and include_cassette) 
            or (chr_type in ('chloroplast', 'mitochondrial', 'other') and include_other)):
            chosen_chromosomes[chr_type].append(chromosome)
    if output_dict:     return dict(chosen_chromosomes)
    else:               return sum(chosen_chromosomes.values(), [])


def _get_plotline_pos(middle_pos, total_width, total_N, N):
    """ Divide a total_width centered at middle_pos into total_N sections, and return the left and right edge of the Nth section. """
    width_per_plot = total_width/total_N
    left_pos = middle_pos - total_width/2 + width_per_plot*N
    right_pos = middle_pos - total_width/2 + width_per_plot*(N+1)
    return left_pos, right_pos


def mutant_positions_and_data(mutant_datasets=[], density=True, colors=None, names='mutants', strands='both', 
                          other_datasets=[], other_density=False, other_colors=None, other_names='other', 
                          bin_size=20000, chromosome_lengths=None, interpolate=False, condense_colorbars=True, 
                          include_scaffolds=False, include_cassette=False, include_other=False):
    """ Plot multiple mutant datasets or other data as position lines or density heatmaps across chromosomes. 

    Mutant_datasets must be a list of mutant_analysis_classes.Insertional_mutant_pool_dataset instances (or a single instance).
    Other_datasets must be a a list of chromosome:position_list dicts (with positions as integers) (or a single dict). 

    Density, colors, names, strands must all be lists of the same length as mutant_datasets, or single values
     (which will be converted to lists filled with those values). 
    The same applies to other_density, other_colors, other_names with respect to other_datasets; 
     otherwise the other_* arguments are the same as the first set, unless otherwise noted. 
    Strands is only needed for mutant_datasets - each mutant dataset will be converted to a chromosome:position_list dict, 
     taking all mutants if the corresponding strands value is 'both', or only +strand/-strand mutants if it's '+' or '-'.

    For each dataset (mutant or other), if the corresponding density value is False, 
      each position will be plotted as a horizontal line (using the corresponding color value, or black if None); 
     if True, the positions will be converted to a density map by binning mutants into ranges (with bin_size bases per range), 
      and plotted as a heatmap (using the corresponding color value as the colormap name, or the default colormap if None).

    The names/other_names values are used for the legend and colorbar labels.

    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default file will be used).

    Interpolate governs whether color blocks are drawn as rectangles or interpolated on density heatmaps. 
    If condense_colorbars, if there is more than one colorbar to be drawn, they'll be condensed into two rows with reduced spacing.

    The include_* True/False arguments govern which chromosome types (in addition to the 17 nuclear chromosomes) should be included.
    """
    ### parse the arguments, set defaults, etc
    # for the arguments that should be lists, if a single value is given, change to a correct-length list filled with that value
    if isinstance(mutant_datasets, mutant_analysis_classes.Insertional_mutant_pool_dataset):    mutant_datasets = [mutant_datasets]
    if isinstance(other_datasets, dict):                                                        other_datasets = [other_datasets]
    if density in (True,False):                                 density = [density for _ in mutant_datasets]
    if other_density in (True,False):                           other_density = [other_density for _ in other_datasets]
    if colors is None or isinstance(colors,str):                colors = [colors for _ in mutant_datasets]
    if other_colors is None or isinstance(other_colors,str):    other_colors = [other_colors for _ in other_datasets]
    if isinstance(names,str):                                   names = [names for _ in mutant_datasets]
    if isinstance(other_names,str):                             other_names = [other_names for _ in other_datasets]
    if isinstance(strands,str):                                 strands = [strands for _ in mutant_datasets]
    imshow_kwargs = {} if interpolate else {'interpolation': 'none'}
    # set sensible default values for None colors (depends on whether it's color or colormap)
    old_colors, old_other_colors = colors, other_colors
    colors, other_colors = [], []
    for old_color_list,new_color_list,density_list in [(old_colors,colors,density), (old_other_colors,other_colors,other_density)]:
        for curr_color,curr_density in zip(old_color_list,density_list):
            if curr_color is not None:  new_color_list.append(curr_color)
            else:   # see http://matplotlib.org/examples/pylab_examples/show_colormaps.html for all colormaps; *_r means reverse
                if curr_density:        new_color_list.append('gist_earth')
                else:                   new_color_list.append('black')

    ### get the list of chromosomes to plot
    # get the chromosome lengths from a file, if filename or None was given
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = _get_chromosome_lengths(chromosome_lengths)
    # grab only the chromosome types we want
    all_chromosomes = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, False)
    # sort chromosomes properly by type and name
    all_chromosomes.sort(key=seq_basic_utilities.chromosome_sort_key)
    # MAYBE-TODO for short chromosomes like scaffolds/chlor/mito/cassette, plot in multiple rows? Maybe change the scale to stretch them a bit, since some will probably be smaller than bin_size?  But then I'd have to normalize them separately or something...

    ### get a list of all datasets provided, converted to the same format, along with their extra args
    #  (for mutant_datasets, need to use _get_mutant_positions_from_dataset to convert each dataset 
    #   to a chromosome:position_list dict; other_datasets already come in that form)
    all_datasets_density_color_name = []
    for dataset,d,c,n,s in zip(mutant_datasets, density, colors, names, strands):
        all_datasets_density_color_name.append((_get_mutant_positions_from_dataset(dataset, strand=s), d, c, n))
    for dataset,d,c,n in zip(other_datasets, other_density, other_colors, other_names):
        all_datasets_density_color_name.append((other_dataset, other_density, other_color))
    total_plotline_width = 0.6

    ### plot all the datasets
    all_heatmaps = []
    for N,(curr_dataset,curr_density,curr_color,curr_name) in enumerate(all_datasets_density_color_name):
        # if we're not doing a density plot, just plot the positions as lines
        if not curr_density:
            for chr_N,chromosome in enumerate(all_chromosomes):
                chromosome_length = chromosome_lengths[chromosome]
                left_pos, right_pos = _get_plotline_pos(chr_N, total_plotline_width, len(all_datasets_density_color_name), N)
                position_list = curr_dataset[chromosome]
                mplt.vlines(chr_N, left_pos, chromosome_length)
                # inverting positions, so that position 0 is at the top, not the bottom - easier to read that way
                inverted_positions = [chromosome_length - pos for pos in position_list]
                # only give a label to one chromosome per dataset, so only one shows up in the legend
                mplt.hlines(inverted_positions, left_pos, right_pos, color=curr_color, 
                            label = ('__nolegend__' if chr_N!=0 else curr_name))
        # if we're doing a density plot, bin the positions into counts, and plot as a heatmap.
        else:
            chromosome_bin_count_lists = {}
            for chr_N,chromosome in enumerate(all_chromosomes):
                chromosome_length = chromosome_lengths[chromosome]
                position_list = curr_dataset[chromosome]
                # divide the chromosome into bin_size-sized ranges; there'll be a smaller-than-bin_size chunk left over, 
                #  so add an extra bin for that if it's at least half a bin_size, otherwise make the last bin bigger to compensate.
                bin_edges = [x-.5 for x in range(1, chromosome_length+1, bin_size)]
                last_bin_edge = chromosome_length-.5
                if (last_bin_edge - bin_edges[-1])/bin_size < .5:   bin_edges[-1] = last_bin_edge
                else:                                               bin_edges.append(last_bin_edge)
                bin_count_list, _  = numpy.histogram(position_list, bin_edges)
                # for the last bin in each chromosome, the value should be scaled by the smaller/bigger bin-size...
                bin_count_list[-1] = bin_count_list[-1] / ((bin_edges[-1] - bin_edges[-2])/bin_size)
                chromosome_bin_count_lists[chromosome] = bin_count_list
            # now calculate the overall maximum bin-count over all chromosomes, so that they're all normalized to the same value
            max_count = max([max(chromosome_bin_count_lists[c]) for c in all_chromosomes])
            for chr_N,chromosome in enumerate(all_chromosomes):
                chromosome_length = chromosome_lengths[chromosome]
                left_pos, right_pos = _get_plotline_pos(chr_N, total_plotline_width, len(all_datasets_density_color_name), N)
                # actually draw the heatmap image - this is tricky! Matplotlib tends to assume there's only a single heatmap.
                #  - the reshape call is to make the densities from a 1D array into a 2D Nx1 array, for imshow to work properly
                #  - aspect='auto' is to prevent matplotlib from reshaping the entire plot to fit the image shape
                #  - norm=None, vmin=0, vmax=1 are to prevent matplotlib from normalizing the values to a 0-1 range!
                #    (which would be BAD, because then each chromosome would be normalized separately, and have a different scale)
                #  - DON'T give those a label, so that mplt.legend() will only work on line-plots - colorbars will be separate
                im = mplt.imshow(chromosome_bin_count_lists[chromosome].reshape(-1,1), 
                                 extent=(left_pos,right_pos,0,chromosome_length), cmap=mplt.get_cmap(curr_color), 
                                 aspect='auto', norm=None, vmin=0, vmax=max_count, **imshow_kwargs)
                # save the info to plot colorbars later
                if chr_N==0:    all_heatmaps.append((im, curr_name))

    ### set plot limits, ticks, labels, etc
    # mplt.imshow has an annoying tendency to reset the plot limits to match a single image, so set them sensibly by hand
    # we want the lower y limit to be slightly below 0, so that the first bin doesn't hide behind the axis line (checked by hand)
    edge_space = 1-total_plotline_width
    mplt.xlim(0 - total_plotline_width/2 - edge_space, len(all_chromosomes)-1 + total_plotline_width/2 + edge_space)
    mplt.ylim(-max(chromosome_lengths.values())*.0015, max(chromosome_lengths.values())*1.01)
    # add xticks with chromosome names as labels!  
    #  - if we only have normal chromosomes, display just the numbers.
    #  - if we have different kinds, need to display the full names - rotate them to keep them from overlapping.
    # MAYBE-TODO add the total number of mutants in each chromosome, and the #mutants/length
    mplt.xlabel('chromosome')
    if not any([include_scaffolds, include_cassette, include_other]):
        mplt.xticks(range(len(all_chromosomes)), [x.split('_')[1] for x in all_chromosomes])
    else:
        mplt.xticks(range(len(all_chromosomes)), all_chromosomes, rotation=90)
    # MAYBE-TODO it'd be nice to get rid of the actual ticks and just keep the labels, the ticks only obscure the data
    # modify yticks to be in Mb - and for some reason this resets ylim, so save it and set it back to previous value
    mplt.ylabel('chromosome position (in MB) (chromosomes start at the top)')
    ylim = mplt.ylim()
    yticks = mplt.yticks()[0]
    mplt.yticks(yticks, [x/1000000 for x in yticks])
    mplt.ylim(ylim)
    # TODO the y scale/direction is inconsistent!  If chromosomes start at the top, they should be lined up on top instead of bottom... Ask Martin what the best way is.
    plotting_utilities.remove_half_frame()

    ### add legends
    # normal legend for the line-plots (only if there's more than one dataset - if there's just one, presumably the title says it)
    if len(all_datasets_density_color_name)>1:
        mplt.legend(prop=FontProperties(size='small'))
        # MAYBE-TODO nicer legend formatting?  Get rid of frame, thicker line, line up with colorbars or something?
    # if desired, figure out sensible positioning to put smaller colorbars in two rows rather than have the default big ones
    if condense_colorbars:
        # column sizes
        N_columns = numpy.ceil(len(all_heatmaps)/2)
        col_width = 0.1
        col_padding = 0.02
        # bbox is the main figure shape - has attributes like xmin,xmax,ymin,ymax,height,width
        ax = mplt.gca()
        bbox = ax.get_position()
        # set new axes position: decrease width to make room for colorbars, keep left/bottom/height as before (from bbox)
        ax.set_position([bbox.xmin, bbox.ymin, bbox.width - (col_padding+col_width)*(N_columns-1), bbox.height])
        mplt.draw()
        bbox = ax.get_position()
        # row positions/sizes
        row_padding = 0.1*bbox.height
        row_height = (bbox.height - row_padding) / 2
        row_bottoms = [bbox.ymin + bbox.height/2 + row_padding/2, bbox.ymin]
    # add colorbars for each heatmap, labeled with the dataset name
    for N, (plot, name) in enumerate(all_heatmaps):
        colorbar_kwargs = {}
        # optionally condense the colorbars: put in two rows, reduce spacing
        if condense_colorbars:
            row, column = N % 2, N // 2
            col1_offset = 1 - N_columns * .05
            # add_axes arguments are left,bottom,width,height
            cax = mplt.gcf().add_axes((bbox.xmax + col_padding + col_width*column, row_bottoms[row], 0.012, row_height))
            colorbar_kwargs['cax'] = cax
        c = mplt.colorbar(plot, **colorbar_kwargs)
        c.set_label("%s (per %s)"%(name,seq_basic_utilities.format_base_distance(bin_size)))
        # MAYBE-TODO put the heatmaps on the plot instead of on the side?  There's room...
        # we don't need as many colorbar ticks as matplotlib likes - leave only the full-number ones
        ticks_and_labels = []
        for x in c._ticker()[1]:
            x = float(x)
            if x==int(x):   ticks_and_labels.append( (x, str(int(x)) ) )
        c.set_ticks(zip(*ticks_and_labels)[0])
        c.set_ticklabels(zip(*ticks_and_labels)[1])
        mplt.draw()


def chromosome_density_scatterplot(mutant_dataset, include_scaffolds=True, include_cassette=True, include_other=True, 
                                        chromosome_lengths=None, chloroplast_multiplier=1, mito_multiplier=1):
    """ Make a chromosome length vs mutant number scatterplot. 

    The include_* True/False arguments govern which chromosome types (in addition to the 17 nuclear chromosomes) should be included.

    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default file will be used).

    The *_multiplier arguments can be set to N>1 to make the plot reflect that fact that there are multiple copies 
     of the chloroplast/mitochondrial genomes in the cell, so the "functional" length is N times higher.
    """

    # get the chromosome lengths from a file, if filename or None was given; grab only the chromosome types we want;
    #  apply chloro/mito multipliers
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = _get_chromosome_lengths(chromosome_lengths)
        try:                chromosome_lengths['chloroplast'] *= chloroplast_multiplier
        except KeyError:    pass
        try:                chromosome_lengths['mitochondrial'] *= mito_multiplier
        except KeyError:    pass
    chromosomes_by_type = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, 
                                                   output_dict=True)

    for chr_type in sorted(chromosomes_by_type, key=seq_basic_utilities.chromosome_sort_key):
        chr_list = chromosomes_by_type[chr_type]
        mutants_in_chromosomes = [mutant_dataset.summary.mutants_in_chromosome(c) for c in chr_list]
        if sum(mutants_in_chromosomes):
            mplt.scatter([chromosome_lengths[c] for c in chr_list], mutants_in_chromosomes, 
                         c=[seq_basic_utilities.chromosome_color(c) for c in chr_list], label=chr_type)

    mplt.ylabel('number of mutants in chromosome')
    mplt.xlabel('chromosome length (in Mb)')
    mplt.xticks(mplt.xticks()[0], [x/1000000 for x in mplt.xticks()[0]])
    mplt.legend(loc='lower right', prop=FontProperties(size='medium'))


def chromosome_density_barchart(mutant_dataset, include_scaffolds=True, include_cassette=True, include_other=True, 
                                        chromosome_lengths=None, chloroplast_multiplier=1, mito_multiplier=1):
    """ Make a simple bar-chart of chromosome mutant densities (per kb).

    See chromosome_density_scatterplot docstring for all the arguments - they're the same.
    """

    # get the chromosome lengths from a file, if filename or None was given; grab only the chromosome types we want; 
    #  apply chloro/mito multipliers; sort
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = _get_chromosome_lengths(chromosome_lengths)
        try:                chromosome_lengths['chloroplast'] *= chloroplast_multiplier
        except KeyError:    pass
        try:                chromosome_lengths['mitochondrial'] *= mito_multiplier
        except KeyError:    pass
    all_chromosomes = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, False)
    all_chromosomes.sort(key=seq_basic_utilities.chromosome_sort_key)

    # calculate and plot mutant densities, colored by chromosome type
    mutant_densities = [mutant_dataset.summary.mutants_in_chromosome(c)/chromosome_lengths[c]*1000 for c in all_chromosomes]

    mplt.bar(range(len(all_chromosomes)), mutant_densities, 
             color=[seq_basic_utilities.chromosome_color(c) for c in all_chromosomes], align='center')

    mplt.ylabel('density of mutants in chromosome (per kb)')
    mplt.xticks(range(len(all_chromosomes)), all_chromosomes, rotation=90)
    mplt.xlim(-1, len(all_chromosomes))
    mplt.legend(loc='lower right', prop=FontProperties(size='medium'))


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


def adjacent_distance_histogram(dataset, incl_same_strand=True, incl_opposite_both=True, incl_opposite_separate=True, 
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


def adjacent_dist1_barchart(dataset, incl_same_strand=True, incl_opposite_both=True, incl_opposite_separate=True, 
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


def adjacent_readcount_ratio_plot(dataset, distance_cutoffs, distance_linestyles=None, incl_same_strand=True, 
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
