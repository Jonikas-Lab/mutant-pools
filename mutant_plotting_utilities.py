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
import random
# other packages
import numpy
import matplotlib.pyplot as mplt
from matplotlib.font_manager import FontProperties
# my modules
import general_utilities
import basic_seq_utilities
import plotting_utilities
import mutant_analysis_classes
from mutant_utilities import get_histogram_data_from_positions, get_mutant_positions_from_dataset, get_chromosome_lengths


########################################## Plotting mutant positions #############################################

######### Plotting mutant/other positions over chromosomes as heatmaps/lines (and help functions for that)

def _get_chromosomes_by_type(all_chromosomes, include_scaffolds, include_cassette, include_other, output_dict=False):
    """ Filter all_chromosomes based on the include_* arguments; return a list of only ones of the desired types. 

    If output_dict is True, output a type:chromosome_list dictionary instead of a single chromosome_list.
    """
    chosen_chromosomes = defaultdict(list)
    for chromosome in all_chromosomes:
        chr_type = basic_seq_utilities.chromosome_type(chromosome)
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


def mutant_positions_and_data(datasets=[], dataset_formats=0, density_plots=True, colors=None, names=None, maxes_per_base=None, 
                              strands=None, title='', bin_size=DEFAULT_BIN_SIZE, chromosome_lengths=None, 
                              interpolate=False, condense_colorbars=True, total_plotline_width=0.6, 
                              include_scaffolds=False, include_cassette=False, include_other=False):
    """ Plot multiple datasets (mutant,position,density) across chromosomes, as position lines or density heatmaps.

    Each element in datasets must match the corresponding value in dataset_formats (or the value if only one is given:
     - format 0 - a mutant_analysis_classes.Insertional_mutant_pool_dataset instance
     - format 1 - a chromosome:position_list dict (with positions as integers)
     - format 2 - a chromosome:bin_density_list dict (with bin densities as numbers - generated with bin sizes matching bin_size.)
    If datasets is a single instance of any of these things instead of a list, it'll be treated as a length-1 list.

    Dataset_formats, density_plots, colors, names, and strands must all be lists of the same length as datasets, or single values
     (which will be converted to lists filled with those values):
     * for each dataset, if the corresponding density_plots value is False, 
        each position will be plotted as a horizontal line (using the corresponding color value, or black if None); 
         (this is impossible if format==2 for this dataset - in that case an error will be returned)
       if True, the positions will be converted to a density map by binning mutants into ranges (with bin_size bases per range), 
        and plotted as a heatmap (using the corresponding color value as the colormap name, or the default colormap if None).
     * the maxes_per_base values can be None or numbers, and are only relevant for density heatmap plots:
        - if None, the density heatmaps will be scaled based on the highest value in the dataset 
            (makes sense for values you wouldn't put on a 0-100% scale, like mutant positions, gene positions, etc)
        - if a number X, the density heatmaps will be scaled so that X positions per base is 100% 
            (makes sense for values with a defined 0-100% range, like mappability or GC content)
     * strands is only needed for mutant datasets - each mutant dataset will be converted to a chromosome:position_list dict, 
        taking all mutants if the strands value is None, or only +strand/-strand/opposite-tandem mutants if it's '+'/'-'/'both'.
     * The names values are used as labels for the legend (if density=False) or colorbar (if density=True); 
        if '' is given, that dataset won't be included in the legend or won't have a colorbar
         (use for cases with multiple copies of similar datasets - they should use the same color as a labeled one!)

    The include_* True/False arguments govern which chromosome types (in addition to the 17 nuclear chromosomes) should be included.
    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default file will be used).
    Interpolate governs whether color blocks are drawn as rectangles or interpolated on density heatmaps. 
    If condense_colorbars, if there is more than one colorbar to be drawn, they'll be condensed into two rows with reduced spacing.
    Title will be used for the plot title.
    Total_plotline_width gives the fraction of x-axis space given to the chromosome plots vs spaces between chromosomes.
    """
    ### parse the arguments, set defaults, etc
    # for the arguments that should be lists, if a single value is given, change to a correct-length list filled with that value
    if isinstance(datasets, (dict, mutant_analysis_classes.Insertional_mutant_pool_dataset)):    datasets = [datasets]
    if isinstance(dataset_formats, int):                                   dataset_formats = [dataset_formats for _ in datasets]
    if density_plots in (True,False):                                      density_plots = [density_plots for _ in datasets]
    if colors is None or isinstance(colors,str):                           colors = [colors for _ in datasets]
    if names is None or isinstance(names,str):                             names = [names for _ in datasets]
    if maxes_per_base is None or isinstance(maxes_per_base,(int,float)):   maxes_per_base = [maxes_per_base for _ in datasets]
    if isinstance(strands,str) or strands is None:                         strands = [strands for _ in datasets]
    imshow_kwargs = {} if interpolate else {'interpolation': 'none'}
    # set sensible default values for None colors (different for mutant and other datasets)
    #  - if density is True, use a colorbar (see ~/computers_and_programming/colormaps_all_matplotlib.png; *_r means reverse)
    #  - otherwise use a single color name
    for N,(color,density_plot,d_format) in enumerate(zip(list(colors),density_plots,dataset_formats)):
        if color is None:
            if density_plot:   colors[N]='gist_earth'     
            else:              colors[N]='black'
    # set sensible default values for None names ('mutants' for mutants, more general for other data, since we can't tell)
    for N,(name,d_format) in enumerate(zip(list(names),dataset_formats)):
        if name is None:
            if d_format==0:    names[N]='mutants'
            else:              names[N]='positions'

    ### get the list of chromosomes to plot
    # get the chromosome lengths from a file, if filename or None was given
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
    # grab only the chromosome types we want
    all_chromosomes = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, False)
    # sort chromosomes properly by type and name
    all_chromosomes.sort(key=basic_seq_utilities.chromosome_sort_key)
    # MAYBE-TODO for short chromosomes like scaffolds/chlor/mito/cassette, plot in multiple rows? Maybe change the scale to stretch them a bit, since some will probably be smaller than bin_size?  But then I'd have to normalize them separately or something...

    ### get a list of all datasets provided, converted to the appropriate format for plotting, along with their extra args:
    #    - chrom:position_list format if density_plot is False
    #    - chrom:bin_density_list format if density_plot is True
    all_datasets_and_info = []
    for dataset,d_format,density_plot,color,name,max_per_base,strand in zip(datasets, dataset_formats, density_plots, 
                                                                            colors, names, maxes_per_base, strands):
        # for mutant datasets, use get_mutant_positions_from_dataset to convert each dataset to a chromosome:position_list dict
        if d_format==0:
            dataset = get_mutant_positions_from_dataset(dataset, strand=strand)
            d_format = 1
        # now convert each dataset to the appropriate format based on density_plot and format, or raise exception if impossible
        if (density_plot is False and d_format==1) or (density_plot is True and d_format==2):
            pass
        elif density_plot is False and d_format==2:
            raise ValueError("Dataset %s was provided in density_plot format - cannot plot as lines!"%name)
        elif density_plot is True and d_format==1:
            dataset = get_histogram_data_from_positions(dataset, bin_size, chromosome_lengths, all_chromosomes, 
                                                        special_last_bin=True, merge_last_bin_cutoff=0.5, normalize_last_bin=True)
        all_datasets_and_info.append((dataset, density_plot, color, name, max_per_base))

    ### plot all the datasets
    all_heatmaps = []
    for N,(dataset,density_plot,color,name,max_per_base) in enumerate(all_datasets_and_info):
        # for heatmaps, decide what to use as the maximum value to scale the colorbar: 
        #  if max_per_base was provided, use that, otherwise use the actual maximum bin-count
        #  (need to use maximum bin-count over all chromosomes, so that they're all normalized to the same value!)
        if density_plot:
            if max_per_base is None:
                max_count = max([max(dataset[c]) for c in all_chromosomes if c in dataset])
            else:
                max_count = max_per_base*bin_size

        # plot data for each chromosome
        for chr_N,chromosome in enumerate(all_chromosomes):
            chromosome_length = chromosome_lengths[chromosome]
            left_pos, right_pos = _get_plotline_pos(chr_N, total_plotline_width, len(all_datasets_and_info), N)
            # if we're not doing a density plot, just plot the positions as lines
            if not density_plot:
                mplt.vlines(chr_N, left_pos, chromosome_length)
                position_list = dataset[chromosome]
                # only plot the lines if there are any mutants
                if position_list:
                    # only give a label to one chromosome per dataset, so only one shows up in the legend
                    mplt.hlines(position_list, left_pos, right_pos, color=color, 
                                label = ('__nolegend__' if (chr_N!=0 or name=='') else name))
            # if we're doing a density plot, bin the positions into counts, and plot as a heatmap.
            else:
                # inverting the bin value list so that 0 is at the bottom, since otherwise the scale makes no sense;
                #  if no data for this chromosome, no need to plot anything.
                try:                inverted_data = numpy.array(list(reversed(dataset[chromosome])))
                except KeyError:    continue
                # actually draw the heatmap image - this is tricky! Matplotlib tends to assume there's only a single heatmap.
                #  - the reshape call is to make the densities from a 1D array into a 2D Nx1 array, for imshow to work properly
                #  - aspect='auto' is to prevent matplotlib from reshaping the entire plot to fit the image shape
                #  - norm=None, vmin=0, vmax=1 are to prevent matplotlib from normalizing the values to a 0-1 range!
                #    (which would be BAD, because then each chromosome would be normalized separately, and have a different scale)
                #  - DON'T give those a label, so that mplt.legend() will only work on line-plots - colorbars will be separate
                im = mplt.imshow(inverted_data.reshape(-1,1), 
                                 extent=(left_pos,right_pos,0,chromosome_length), cmap=mplt.get_cmap(color), 
                                 aspect='auto', norm=None, vmin=0, vmax=max_count, **imshow_kwargs)
                # save info image to add colorbars: image, name, and if_normalized (True if max_per_base was given) 
                if chr_N==0 and name!='':  all_heatmaps.append((im, name, (max_per_base is not None)))

    ### set plot limits, ticks, labels, etc
    # it's important to do all this BEFORE colorbars, otherwise it'll apply to the colorbar instead of main axes!
    if title:  mplt.title(title)
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
    mplt.ylabel('chromosome position (in Mb) (chromosomes start at the bottom)')
    ylim = mplt.ylim()
    yticks = mplt.yticks()[0]
    mplt.yticks(yticks, [x/1000000 for x in yticks])
    mplt.ylim(ylim)
    plotting_utilities.remove_half_frame()

    ### add legends
    # normal legend for the line-plots (only if there's more than one dataset - if there's just one, presumably the title says it)
    if len(all_datasets_and_info)>1:
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
    for N, (plot, name, if_normalized) in enumerate(all_heatmaps):
        colorbar_kwargs = {}
        # optionally condense the colorbars: put in two rows, reduce spacing
        if condense_colorbars:
            row, column = N % 2, N // 2
            col1_offset = 1 - N_columns * .05
            # add_axes arguments are left,bottom,width,height
            cax = mplt.gcf().add_axes((bbox.xmax + col_padding + col_width*column, row_bottoms[row], 0.012, row_height))
            colorbar_kwargs['cax'] = cax
        c = mplt.colorbar(plot, **colorbar_kwargs)
        c.set_label("%s"%name + ("" if if_normalized else " (per %s)"%(basic_seq_utilities.format_base_distance(bin_size))))
        # MAYBE-TODO put the heatmaps on the plot instead of on the side?  There's room...
        # if if_normalized, scale the ticks to 0-100%
        if if_normalized:
            # c.vmin and c.vmax are the low/high of the colorbar
            assert c.vmin==0, "Colorbar start isn't at 0, can't normalize to 0-100% properly!"
            ticks_and_labels = [(c.vmax*fraction, "%d%%"%(fraction*100)) for fraction in (0, 0.25, 0.5, 0.75, 1)]
        # otherwise, we don't need as many colorbar ticks as matplotlib likes - leave only the full-number ones
        else:
            ticks_and_labels = []
            # c._ticker()[1] gives the tick positions
            for x in c._ticker()[1]:
                x = float(x)
                if x==int(x):   ticks_and_labels.append( (x, str(int(x)) ) )
        c.set_ticks(zip(*ticks_and_labels)[0])
        c.set_ticklabels(zip(*ticks_and_labels)[1])
        mplt.draw()
    # TODO it'd be nice if this actually returned the main axes object... (and maybe a list of colorbar axes too?)  And the numpy.histogram results for all the datasets might be nice too!


######### Hotspot/coldspot related plots

def _lowest_cutoff_matched(val, cutoff_list):
    """ Given x and a list of cutoffs, return how many of them x is lower than. """
    return sum(val<cutoff for cutoff in cutoff_list)


def get_hot_cold_spot_colors(pvalue_data, side_data, window_size, window_offset, pval_cutoffs=[0.05, 0.01, 0.005, 0.001], 
                             min_zero=False, print_info=True):
    """ Transform p-value and side (low/high) chrom:value_list dicts into plottable chrom:value lists.

    Both arguments should be chromosome:list_of_window_values, with any number of windows, with given window size and offset. 
      - the p-values should be 0-1, and you should probably get FDR correction done on them first; 
      - the sides should be -1 for coldspots and 1 for hotspots. 
    These can be generated by mutant_simulations.find_hot_cold_spots (return values 1 and 3)

    Pval_cutoffs should be a list of 0-1 values of any length, sorted reverse.
    The output will be a chrom:heatmap_value_list, with the value 0 for p-values over pval_cutoffs[1], 
     1 or -1 for p-values under pval_cutoffs[1] but over [2], and so on, with the maximum being len(pval_cutoffs), 
     and the minimum the negative of that.  Negative values indicate coldspots, and positive hotspots. 
    If min_zero is True, len(pval_cutoffs) will be added to each value, to make the lowest possible value 0
     (can be easier for plotting).
    Suitable for plotting with mutant_positions_and_data, format 2, using a blue-white-red colormap or similar, 
     IF the offset and the remainder at the end of the chromosome are relatively small!  Otherwise may end up messy.

    If print_info, print a line giving the number of somewhat significant results, just to give a general idea.
    """
    # TODO modify mutant_positions_and_data to deal with negative values, and to customize colorbar ticklabels for p-values! Add that to docstring both there and here when done
    # TODO implement refactoring the whole list to deal with offsets so it gets plotted correctly, if we see anything worth plotting (right now the positions for different offsets are the same, and the missing chunks at chromosome starts/ends are ignored)
    format_bp = lambda x: basic_seq_utilities.format_base_distance(x, False)    # the False is to not approximate
    color_value_dict = {}
    N_significant = 0
    for chrom, pvalues in pvalue_data.items():
        sides = side_data[chrom]
        color_values = []
        for N,(pvalue,side) in enumerate(zip(pvalues,sides)):
            curr_val = _lowest_cutoff_matched(pvalue, pval_cutoffs)
            if curr_val>0 and print_info:
                N_significant += 1
            curr_val *= side
            if min_zero:    curr_val += len(pval_cutoffs)
            color_values.append(curr_val)
        color_value_dict[chrom] = color_values
    if print_info:
        print "%s results below adjusted p-value %s for data with %s window (offset %s)"%(N_significant, pval_cutoffs[0], 
                                                                               format_bp(window_size), format_bp(window_offset))
    return color_value_dict
    # TODO should unit-test this!


def plot_hot_cold_spot_hlines(hc_spot_list, pval_cutoffs, all_chromosomes=None, min_offset=0, max_offset=0.4, 
                              linewidth=2, min_height=10000, pval_cutoff_colors=None, chrom_numbers=None):
    """ Plot hot/cold spots as horizontal lines, colored by pvalue by chromosome position. 
    
    Hc_spot_list should be a list of (chrom, start, end, pvalue, kind, window_offset) tuples 
     as generated by mutant_simulations.get_hot_cold_spot_list (with kind being 'hotspot' or 'coldspot', 
     and window_offset along with start-end size used to separate potentially overlapping data into lines),
    Pval_cutoffs should be a reverse-sorted list of cutoffs, 0-1, to be used when coloring the data. 

    Min/max offset give the x-axis range relative to the middle of each chromosome over which to spread all the data
     (1 is the distance between two chromosomes on the x axis, so the full range would be about -.4 to .4)
    Min_height is the minimum height 

    Pval_cutoff_colors should be a kind:color_list dict, with each color_list the length of pval_cutoffs, 
     or default if None.

    All_chromosomes can be a list of all chromosomes to include on the plot, 
     in case you want to leave space for any that DON'T have any hc_spot_list entries.
    Chrom_numbers should be a chromosome_name:number dictionary to give each chromosome a unique position on the x axis; 
     if None, we'll just sort them all sensibly and assign numbers that way. 
    """
    if pval_cutoff_colors is None:
        if len(pval_cutoffs) == 3:
            pval_cutoff_colors = {'hotspot': 'salmon red darkred'.split(), 'coldspot': 'skyblue steelblue darkblue'.split()}
        else:
            raise Exception("If pval_cutoffs isn't length 3, must provide pval_cutoff_colors, no default!")
    if all_chromosomes is None:
        all_chromosomes = list(set(chrom for (chrom, start, end, pvalue, kind, offset) in hc_spot_list))
    else:
        all_chromosomes = list(all_chromosomes)
    if chrom_numbers is None:
        all_chromosomes.sort(key=basic_seq_utilities.chromosome_sort_key)
        print "Sorted chromosomes: %s"%(', '.join(all_chromosomes))
        chrom_numbers = {chrom:N for (N,chrom) in enumerate(all_chromosomes)}
    # figure out how many potentially overlapping (window_size, window_offset) sets we have, 
    #  give each of them a separate x axis position
    all_lines = sorted(set((end-start, offset) for (chrom, start, end, pvalue, kind, offset) in hc_spot_list))
    x_offset_range = max_offset-min_offset
    x_offset_per_line = x_offset_range/len(all_lines)
    x_line_offsets = {line: (min_offset)+x_offset_per_line*N for (N,line) in enumerate(all_lines)}
    for (chrom, start, end, pvalue, kind, offset) in hc_spot_list:
        # plot only the lines with a pvalue that matches at least the highest cutoff
        cutoff_matched = _lowest_cutoff_matched(pvalue, pval_cutoffs)
        if cutoff_matched>0:
            x_offset = x_line_offsets[(end-start,offset)]
            # if min_height is specified, make sure the line height is at least that
            if min_height is not None:
                height = end-start
                missing_height = min_height - height
                if missing_height>0:
                    start, end = start - missing_height/2, end + missing_height/2
            # line color depends on lowest pvalue cutoff matched
            color = pval_cutoff_colors[kind][cutoff_matched-1]
            mplt.vlines(chrom_numbers[chrom]+x_offset, start, end, color=color, linewidth=linewidth)


######### Gap size plots

def get_gap_sizes(dataset):
    """ Given a dataset, get a list of all the gap sizes between mutant positions, UNSORTED. 
    
    Go over each chromosome, take each pair of adjacent mutants (regardless of strand), get gap size.

    Dataset can be either a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a chrom:position_list dict.
    """
    # convert dataset to chrom:position_list dict if needed
    if isinstance(dataset, mutant_analysis_classes.Insertional_mutant_pool_dataset):
        dataset = get_mutant_positions_from_dataset(dataset, strand=None)
    all_gap_sizes = []
    # go over each chromosome and give the gaps between each adjacent pair of mutants
    for chrom,pos_list in dataset.items():
        pos_list.sort()
        for pos1,pos2 in zip(pos_list, pos_list[1:]):
            all_gap_sizes.append(pos2-pos1)
    return all_gap_sizes


def gap_size_QQ_plot(reference_dataset, simulated_datasets=None, 
                     N_simulations=None, mappable_positions_20bp=None, mappable_positions_21bp=None, 
                     logscale=True, different_colors=False, markersize=6, plot_padding=(0.1, 0.1)):
    """ Quantile-quantile plot of gap sizes between mutants in real vs simulated reference_dataset (repeated N times). 

    Provide either a list of simulated_datasets (chrom:pos_list dicts), OR information to generate them:
     N_simulations, and the two mappable_positions_* args (output of mutant_simulations.genome_mappable_insertion_sites).
    
    """
    reference_gap_sizes = get_gap_sizes(reference_dataset)
    # Problem: since the genome's separated into chromosomes, simulated datasets with the same number of mutants
    #   can have different numbers of gap sizes!  2 mutants on chr1 and chr2 - 2 gaps; 4 and 0 - 3 gaps. 
    #  This difference can't be more than the number of chromosomes, and likely much less. 
    #  So just remove 10 random gap sizes just in case.  MAYBE-TODO come up with a better solution?
    random.shuffle(reference_gap_sizes)
    reference_gap_sizes = sorted(reference_gap_sizes[:-10])
    max_gap = max(reference_gap_sizes)
    # if simulated datasets were provided, just use them; otherwise generate data to do the simulations
    if simulated_datasets is not None:
        N_simulations = len(simulated_datasets)
    else:
        N_mutants = len(reference_dataset)
        fraction_20bp = mutant_simulations.get_20bp_fraction(reference_dataset)
        # TODO mutant_simulations.get_20bp_fraction only works on a full mutant dataset, not on a chrom:pos_list dict - either require reference_dataset to be real, or add an option to give that separately
    if different_colors:    plot_kwargs = {}
    else:                   plot_kwargs = {'color': 'blue'}
    for N in range(N_simulations):
        # if simulated datasets were provided, use the next one;
        #  otherwise make a new one for each repetition (and discard it afterward, to avoid taking up memory)
        if simulated_datasets is not None:
            simulated_dataset = simulated_datasets[N]
        else:
            simulated_dataset = mutant_simulations.simulate_dataset_mappability(N_mutants, fraction_20bp, 
                                                                                mappable_positions_20bp, mappable_positions_21bp)
        simulated_gap_sizes = get_gap_sizes(simulated_dataset)
        # make sure the number of simulated gap sizes matches the number of reference ones
        #   (see "Problem" comment section up top for why this can happen)
        if len(simulated_gap_sizes) < len(reference_gap_sizes):
            raise Exception("Fewer simulated than real gap sizes! Remove more from the real ones?")
        random.shuffle(simulated_gap_sizes)
        simulated_gap_sizes = sorted(simulated_gap_sizes[:len(reference_gap_sizes)])
        max_gap = max(max_gap, max(simulated_gap_sizes))
        if N==0:    label = "real vs simulated,\nrepeated %s times"%N_simulations
        else:       label = '__nolegend__'
        mplt.plot(reference_gap_sizes, simulated_gap_sizes, '.', markeredgecolor='none', label=label, markersize=markersize, **plot_kwargs)
    mplt.title("Real vs simulated gap sizes between insertions (%s scale),"%('log' if logscale else 'linear')
               +"\nquantile-quantile plot (comparing sorted gap-size lists)"
               +"\n(simulation taking mappability into account, repeated %s times)"%N_simulations)
    mplt.xlabel('Real gap sizes')
    mplt.ylabel('Simulated gap sizes')
    if logscale:
        # use symlog instead of log to deal with gap size 0 - the range between -1 and 1 will be linear
        mplt.xscale('symlog',linthreshx=1)
        mplt.yscale('symlog',linthreshy=1)
        # MAYBE-TODO are those sensible edges?  
        #  Could come up with some way of taking plot_padding[0] into account, but it's tricky on a symlog plot
        plot_edges = (-0.5, 10**((1+plot_padding[1])*math.log10(max_gap)) )
    else:
        plot_edges = (-plot_padding[0]*max_gap, (1+plot_padding[1])*max_gap)
    # make the plot square even if the data isn't; plot a diagonal line to show what identical datasets would look like
    mplt.plot(plot_edges, plot_edges, c='grey', label='if both were identical')
    mplt.xlim(plot_edges)
    mplt.ylim(plot_edges)
    mplt.legend(loc=2)
    # TODO those plots should be square!  How do I square an axes instance?


######### Chromosome mutant density plots

def chromosome_density_scatterplot(mutant_dataset, include_scaffolds=True, include_cassette=True, include_other=True, 
                                        chromosome_lengths=None, chloroplast_multiplier=1, mito_multiplier=1):
    """ Make a chromosome length vs mutant number scatterplot. 

    The include_* True/False arguments govern which chromosome types (in addition to the 17 nuclear chromosomes) should be included.

    Chromosome_lengths can be either a chromosome:length dict, or the name of a genome fasta file to extract them from
     (if None, the default file will be used).

    The *_multiplier arguments can be set to N>1 to make the plot reflect that fact that there are multiple copies 
     of the chloroplast/mitochondrial genomes in the cell, so the "functional" length is N times higher.
    """
    # TODO add option to only count effective (mappable) length!

    # get the chromosome lengths from a file, if filename or None was given; grab only the chromosome types we want;
    #  apply chloro/mito multipliers
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
        try:                chromosome_lengths['chloroplast'] *= chloroplast_multiplier
        except KeyError:    pass
        try:                chromosome_lengths['mitochondrial'] *= mito_multiplier
        except KeyError:    pass
    chromosomes_by_type = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, 
                                                   output_dict=True)

    mutants_in_chromosomes_all = {}
    for chr_type in sorted(chromosomes_by_type, key=basic_seq_utilities.chromosome_sort_key):
        chr_list = chromosomes_by_type[chr_type]
        mutants_in_chromosomes = [mutant_dataset.summary.mutants_in_chromosome(c) for c in chr_list]
        mutants_in_chromosomes_all.update(dict(zip(chr_list, mutants_in_chromosomes)))
        if sum(mutants_in_chromosomes):
            mplt.plot([chromosome_lengths[c] for c in chr_list], mutants_in_chromosomes, 'o', 
                      color=basic_seq_utilities.CHROMOSOME_TYPE_COLORS[chr_type], label=chr_type)

    max_length = max(chromosome_lengths.values())
    max_N_mutants = max(mutants_in_chromosomes_all.values())
    mplt.xlabel('chromosome length (in Mb)')
    mplt.xticks(mplt.xticks()[0], [x/1000000 for x in mplt.xticks()[0]])
    mplt.xlim(-max_length*.05, max_length*1.1)
    mplt.ylabel('number of mutants in chromosome')
    mplt.ylim(-max_N_mutants*.05, max_N_mutants*1.1)
    mplt.legend(loc='lower right', prop=FontProperties(size='medium'))


def chromosome_density_barchart(mutant_dataset, include_scaffolds=True, include_cassette=False, include_other=True, 
                                        chromosome_lengths=None, chloroplast_multiplier=1, mito_multiplier=1):
    """ Make a simple bar-chart of chromosome mutant densities (per kb).

    See chromosome_density_scatterplot docstring for all the arguments - they're the same.
    """
    # TODO add option to only count effective (mappable) length!

    # get the chromosome lengths from a file, if filename or None was given; grab only the chromosome types we want; 
    #  apply chloro/mito multipliers; sort
    if chromosome_lengths is None or isinstance(chromosome_lengths, str):
        chromosome_lengths = get_chromosome_lengths(chromosome_lengths)
        try:                chromosome_lengths['chloroplast'] *= chloroplast_multiplier
        except KeyError:    pass
        try:                chromosome_lengths['mitochondrial'] *= mito_multiplier
        except KeyError:    pass
    all_chromosomes = _get_chromosomes_by_type(chromosome_lengths.keys(), include_scaffolds, include_cassette, include_other, False)
    all_chromosomes.sort(key=basic_seq_utilities.chromosome_sort_key)

    # calculate and plot mutant densities, colored by chromosome type
    mutant_densities = [mutant_dataset.summary.mutants_in_chromosome(c)/chromosome_lengths[c]*1000 for c in all_chromosomes]

    mplt.bar(range(len(all_chromosomes)), mutant_densities, 
             color=[basic_seq_utilities.chromosome_color(c) for c in all_chromosomes], align='center', linewidth=0)
    # MAYBE-TODO add error bars based on total mutant number?
    # MAYBE-TODO add a line/something showing the max mutants/20kb value for each chromosome?  But that would require more data processing (similar to what was done in mutant_positions_and_data)

    mplt.ylabel('mutants per kb')
    mplt.xticks(range(len(all_chromosomes)), all_chromosomes, rotation=90)
    mplt.xlim(-1, len(all_chromosomes))


######################################### Plotting gene-related data ###########################################

### number of genes with 1+/2+/etc mutants vs number of mutants (randomly chosen mutant subsets)

def genes_with_N_mutants(dataset, step_size=100, max_N_mutants=3, N_mutants_colors=None, repeat_N_times=100, 
                         total_genes=None, mappable_percent=None, plot_full_count=True, print_repeat_progress=10):
    """ Plot % of all genes with N mutants for different-sized random subsets of dataset.

    Plot % of all genes that have between 1 and max_N_mutants mutants (separate line for each N); 
     N_mutants_colors should be a list of length max_N_mutants, giving the colors to plot the genes with each N mutants.
    Total_genes is the total number of genes in the genome; if None, it'll be extracted from dataset if possible, 
     or else the default value from the Phytozome chlamy v4.3 genome will be used (17114).

    Plot points for mutant subset sizes between 0 and almost the full dataset size, in increments of step_size; 
     the last point won't be the full dataset size, but the closest lower number divisible by step_size.

    Dataset should be a mutant_analysis_classes.Insertional_mutant_pool_dataset instance, 
     or a list of mutants (mutant_analysis_classes.Insertional_mutant instances).

    Redo repeat_N_times, with different random subsets each time, to see the full representation.
    If plot_full_count is True, also plot dots for the full dataset (in case its size is not divisible by step_size), 
     in a different style (black borders around the dots).

    Mappable_percent is the approximate % of mutants that are mappable, just to put in the xlabel; if None, don't print it. 
    """
    # TODO update to show simulated data in parallel to real data, once we have simulated data!  (Simulated data can go to higher total mutant numbers than we have in dataset)  Simulated lines should be paler colors or different markers or something.

    # nice default color schemes only defined for some max_N_mutants values
    if max_N_mutants<=3 and N_mutants_colors is None:
        N_mutants_colors = 'magenta cyan orange'.split()

    # default total_genes: extract from dataset if possible; 
    #  if not (dataset missing that attribute, dataset is a mutant list, or dataset value is None as well), use hardcoded default. 
    if total_genes is None:
        try:                    total_genes = dataset.total_genes_in_genome
        except AttributeError:  pass
        if total_genes is None: total_genes = 17114

    # Plot it all several times with new random mutant subsets, to make sure we have a good coverage of the random space
    for repeat in range(repeat_N_times):
        if print_repeat_progress:
            if not repeat % print_repeat_progress:  
                print "repeat %s/%s..."%(repeat, repeat_N_times)
        gene_counts_by_Nmutants = mutant_simulations.gene_counts_for_mutant_subsets(dataset, max_N_mutants, step_size)
        for N_mutants,gene_counts in gene_counts_by_Nmutants.items():
            # use custom colors if given, otherwise let mplt choose colors
            plot_kwargs = {} if N_mutants_colors is None else {'color': N_mutants_colors[N_mutants-1]}
            # only label the lines in the first repeat, to avoid 100 repeats in the legend!
            if repeat==0:                   plot_kwargs['label'] = "genes with %s+ mutants"%N_mutants
            mplt.plot(gene_counts, '.', linewidth=0, **plot_kwargs)
    # add a separate plot set for the actual full mutant count, with different style (black border)
    if plot_full_count:
        gene_counts_by_Nmutants = mutant_simulations.gene_counts_for_mutant_subsets(dataset, max_N_mutants, 
                                                                                single_subset_size=len(dataset))
        for N_mutants,gene_counts in gene_counts_by_Nmutants.items():
            plot_kwargs = {} if N_mutants_colors is None else {'color': N_mutants_colors[N_mutants-1]}
            mplt.plot(len(dataset)/step_size, gene_counts, '.', linewidth=0, markeredgecolor='black',**plot_kwargs)
    # add a line at the total gene number - TODO this doesn't work right... figure out sensible xmin/xmax values, reset xlim
    mplt.legend(loc='upper left', prop=FontProperties(size='medium'))

    mplt.title('Number of genes hit vs number of mutants sequenced\n(plotted for %s random mutant subsets)'%repeat_N_times)
    mplt.ylabel("Number of genes hit (out of %s total chlamy nuclear genes)"%total_genes)
    mplt.yticks(mplt.yticks()[0], ["%i (%.0f%%)"%(x, x*100/total_genes) for x in mplt.yticks()[0]])
    mplt.xlabel("Number of mutants (randomly chosen out of %s total)\n"%len(dataset) 
                +"(counting only mappable unique genomic mutants%s)"%(' - about %s%%'%mappable_percent if mappable_percent else ''))
    mplt.xticks(mplt.xticks()[0], ["%i"%(x*100) for x in mplt.xticks()[0]])


############################################## Plotting readcounts ##############################################

######### Help functions

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
    # MAYBE-TODO this function may be applicable to more than readcounts - if so, move it up?  Should probably remove the readcounts_only/etc options in that case.  So really that would be a different function...


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
