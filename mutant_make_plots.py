#! /usr/bin/env python
"""
For a list of RNAi library screen deepseq count or growth rate input files, plot the coverages/distributions (a line per input dataset showing the number of shRNAs with given count/growthrate values, on one graph) and/or the correlations (a subplot for each pair of datasets showing the scatterplot correlation between them, optionally also giving Spearman/Pearson correlation coefficients). 
Note that this can be very slow with many/large datasets - use the -t/-s options to generate only the graphs you actually want.
 --Weronika Patena, 2010-2011

USAGE: mutant_make_plots.py [options] infile
"""
# TODO fix docstring!

### basic libraries
import os
from collections import defaultdict
### other packages
from numpy import isnan, isinf, isneginf
from numpy import corrcoef
from scipy.stats.stats import spearmanr
# TODO might want to import some of the bigger libraries only after the option-parsing works to save time...
# switch backend to Agg - this is only necessary on some computers - seems to fix some weird display problem
import matplotlib; matplotlib.use('Agg')    
# then import matplotlib.pyplot normally
import matplotlib.pyplot as mplt

### my modules
import general_utilities
from mutant_analysis_classes import *


############# A LOT OF THIS IS OLD CODE FROM A SIMILAR PROGRAM, TODO CHANGE OR DELETE

def get_sample(dataset, all_mutants, replace_zeros = 0.1):
    """ Return the data list for the sample, in the order of all_mutants, with low and missing values set to 0.1. """
    sample_data = []
    for mutant in all_mutants:
        try:                
            value = dataset[mutant]
            if value < replace_zeros:   value = replace_zeros
        except KeyError:    
            value = replace_zeros
        sample_data.append(value)
    return sample_data

def remove_nan_inf_from_two_samples(sample1, sample2, min_val=None):
    """ Given two float lists, return two lists with any elements that are nan/inf in EITHER list removed from both. """
    assert len(sample1)==len(sample2), "two samples aren't the same length! %s, %s."%(len(sample1),len(sample2))
    new_sample1, new_sample2 = [], []
    for x,y in zip(sample1,sample2):
        if isnan(x) or isnan(y) or isinf(x) or isinf(y) or isneginf(x) or isneginf(y) or x<min_val or y<min_val:
            continue
        new_sample1.append(x)
        new_sample2.append(y)
    return new_sample1, new_sample2

def sort_oligo_list(oligo_list):
    """ Return the reference oligo list, sorted by group (CD45 etc, largest to smallest). """
    # is there a reason for this?  Oh, just to make sure the smallest group gets plotted last so the dots are still visible
    oligo_groups = [[],[],[],[]]
    for oligo_name in oligo_list:
        if oligo_name.count('CD45'):    oligo_groups[0].append(oligo_name)
        elif oligo_name.count('CD19'):  oligo_groups[1].append(oligo_name)
        elif oligo_name.count('TETr'):  oligo_groups[2].append(oligo_name)
        elif oligo_name.count('TetR'):  oligo_groups[2].append(oligo_name)
        else:                           oligo_groups[3].append(oligo_name)
    oligo_groups.sort(key = lambda a: len(a), reverse=True) # sort the groups by size, largest to smallest
    return reduce(lambda a,b:a+b, oligo_groups) # return a concatenation of all the groups into one list

# make color of each dot different for CD45, CD19 and TetR shRNAs:  CD45 - blue, CD19 - green, TetR - red
def find_color(oligo_name):
    if oligo_name.count('CD45'):    return 'b'
    elif oligo_name.count('CD19'):  return 'r'
    elif oligo_name.count('TETr'):  return 'g'
    elif oligo_name.count('TetR'):  return 'g'
    else:                           return 'k'

### Do a standard coverage plot for all the deepseq samples in input
# INPUT: all_data is a list of deepseq_dataset objects (either counts or ratios)
# TODO this currently doesn't work at all!!
def plot_all_distributions_sorted(all_data,plot_scale,figname,value_type):
    library_oligos = all_data[0].data.keys()
    all_samples = {}
    for dataset in all_data:
        sample = get_sample(dataset,library_oligos)
        # the reason the sorting doesn't work is because there are 'nan' values - change them all to 0.
        sample = [x if not isnan(x) else 0 for x in sample]
        sample.sort()
        all_samples[dataset.name] = sample
    sample_names = [dataset.name for dataset in all_data]
    fig = mplt.figure(dpi=200)
    if plot_scale=='log': mplt.semilogy()
    mplt.title('%s PLOT, %s SCALE.\n For each deep-sequencing sample, the %s for each oligo, \nsorted independently in increasing order.'%(os.path.basename(figname.upper()),plot_scale.upper(),value_type), horizontalalignment='center',verticalalignment='bottom')
    mplt.xlabel('oligos in the library, sorted independenty for each sample')
    mplt.ylabel(value_type)
    for name in sample_names:
        mplt.plot(all_samples[name],',',label=name)
    # TODO move the legend somewhere off the actual figure area!  Under should be good.  Some info: http://matplotlib.sourceforge.net/users/legend_guide.html#plotting-guide-legend, http://old.nabble.com/Stopping-Legend-From-Overlapping-the-Graph-td24213554.html, http://matplotlib.sourceforge.net/users/plotting/legend.html#legend-location, http://matplotlib.sourceforge.net/examples/pylab_examples/legend_demo3.html
    #mplt.legend(bbox_transform=mplt.gcf().transFigure,bbox_to_anchor=(0,0,1,1),loc=1)
    #mplt.legend(bbox_transform=mplt.gcf().transFigure,loc=(-1,-1))
    # TODO make the legend font smaller?  This is apparently really hard.
    mplt.legend(loc=2,numpoints=10,handlelength=1.5,handletextpad=0.2)
    return fig

# TODO another type of distribution plot it would be good to have:  a basic histogram-type thing showing the number of shRNAs with each count, listing/showing the mean/median and stdev of each sample, and maybe optionally with a pale line for the matching normal distribution 
def plot_all_distributions_histogram(all_data,plot_scale,figname,value_type):
    pass


def _get_single_clean_sample(all_sample_data, which_sample, min_value=None):
    #print len(all_sample_data), len(all_sample_data[0])
    #print all_sample_data[0]
    #print "grabbing sample %s, with min_value %s"%(which_sample, min_value)
    sample = [float(data_line[which_sample]) for data_line in all_sample_data]
    if verbose:
        print "Parsing sample %s, length %s"%(which_sample,len(sample))
    #print "original:", sample[:10]+sample[-3:]
    isweird = lambda x: bool(isnan(x) or isinf(x) or isneginf(x))
    if min_value is not None:
        sample = [x if (x>min_value or isweird(x)) else min_value for x in sample]
    # TODO for numeric values that just can't be plotted on a given plot (like inf, or 0 for logscale) change them to a min/max value, but also COLOR/SHAPE THEM DIFFERENTLY and put a legend saying that those dots actually represent 0/inf/whatever, and are just plotted as a plottable value.
    # TODO actually for changing/removing zero/negative values from a log plot I could just use nonposx/nonposy arguments to xscale/yscale functions (set to 'mask' to remove or 'clip' to change to small values) instead of doing this by hand...
    #print "cleaned up:", sample[:10]+sample[-3:]
    return sample

def _make_sample_label_text(sample_name, sample, datatype, min_value=0, end_spaces=0, top_newlines=0, end_newlines=0):
    # LATER-TODO should probably adjust the position in some better way than this hack with spaces and dots at the end - maybe just instead of using xlabel/ylabel, add another row/column of subplots and put text in those?
    text = '\n'*top_newlines + sample_name
    if end_spaces:  text+=' '*end_spaces + '.'
    if datatype=='read counts':  
        N_reads = sum([x for x in sample if x>=1])
        text += '\n(%.1fM reads)'%(N_reads/1000000.0)
        if end_spaces:  text+=' '*end_spaces + '.'
    text += '\n(%s mutants)'%len([x for x in sample if x>min_value])
    if end_spaces:  text+=' '*end_spaces + '.'
    text += '\n'*end_newlines
    # TODO is this #mutants meaningful for growth rates, and does it work sensibly?
    return text

### Visualize correlations between all pairs in a set of N deepseq samples 
# Generates a figure with a subplot for each pair, labeled, with correlation coefficient.
# INPUT: all_data is a name->data dictionary, where data is a (details,total_counts) tuple,
#  in which details is a library_oligo_name->deepseq_counts dictionary.
def plot_all_correlations(all_data, plot_scale, figname, print_correlation=False, min_value=0.1, mutants_to_color=None):
    """ _________
    all_data will be an (sample_data, sample_headers, mutant_data) tuple.
    """
    # TODO docstring
    if figname.lower().count('growthrate'):     datatype = 'growth rates'
    elif figname.lower().count('td'):           datatype = 'doubling times'
    elif figname.lower().count('doubling'):     datatype = 'doubling times'
    elif figname.lower().count('count'):        datatype = 'read counts'
    else: datatype = 'unknown'
    # TODO what about TD changes or growth rate changes, do those need their own datatypes?
    sample_data, sample_names, mutant_data = all_data
    N_samples = len(sample_names)
    if mutants_to_color:    colors = [mutants_to_color[tuple(mutant_line[:3])] for mutant_line in mutant_data]
    else:                   colors = ['k' for mutant_line in mutant_data]
    # TODO might want plotsize to be an option, along with default dotsize, and maybe dpi?  Does the dpi even do anything?
    plotsize = 4
    dotsize = 1
    # TODO does the dpi here even do anything??  not sure - as far as the final image it's the savefig dpi that matters.
    # TODO make all the subplots square!!  Using "(plotsize*(N_samples+0.1), plotsize*N_samples)" is close but not exact
    fig = mplt.figure(figsize=(plotsize*(N_samples+0.1), plotsize*N_samples), dpi=300)
    if print_correlation=='none':        correlation_info = ''
    elif print_correlation=='spearman':  correlation_info = '\n\'corr\' = Spearman rank correlation.'
    elif print_correlation=='pearson':   correlation_info = '\n\'corr\' = Pearson\'s correlation coefficient.'
    fig_title = '%s PLOT, %s SCALE. \nFor each pair of samples in the set, a plot comparing %s of each mutant.%s'%(os.path.basename(figname.upper()), plot_scale.upper(), datatype, correlation_info)
    # TODO more explanation?  Also note that more explanatory text could be moved into the empty space.
    # TODO alternatively, the correlation info could be moved into the empty space instead of sticking it on the plots (or put it on the plots, and ALSO put it in the empty space in a table with more info)
    # TODO should I in fact be ignoring mutants not found in one sample? I'd rather plot them...
    # TODO useful idea: make all the dots on the scatterplots have sizes based on # mutants! (optionally). Like, instead of plotting the dots on top of each other, first make a dictionary by dot coordinates, that counts the number of dots with each coordinate pair (there's going to be a lot with (1,1) and (0,1) etc, and presumably just singles with higher numbers), and then plot each coordinate pair just once, with size based on the count!  This should be OPTIONAL.
    #   should probably put a scale somewhere, that gives the size of a dot of 1, 2 or 5 or something, 10, and maybe of the largest dot present on the graph...  Well, on correlation plots I have plenty of room on the lower-left side for scales and such.
    #   will probably need an exception for (0,0), don't really need to plot that
    # TODO if mutants_to_color was used, print a legend listing the colors somewhere in the empty lower-left corner!
    # if it's a logscale plot of readcounts, replace 0's with 0.1
    # TODO what to do if it's a logscale plot of growth rates?  What kinds of growth rate values do we normally see, what's a good minimum?
    # TODO remember to change the scale from 0.1 to 0 too!
    min_value = min_value if plot_scale=='log' else None
    min_val_to_plot = None
    grey_out_min = None
    # MAYBE-TODO min_val_to_plot and grey_out_min could also be a list with a different minimum for each sample...
    print "min_value: %s; min_val_to_plot: %s, grey_out_min: %s"%(min_value, min_val_to_plot, grey_out_min)
    # TODO minimize_tick_labels and same_scale_both_plots should be command-line options
    minimize_tick_labels = False
    same_scale_both_plots = True
    # TODO min_value and min_val_to_plot should probably be merged?  Just one min_value and one True/False variable specifying whether values under min_value should be set to min_value or removed.
    need_title = True
    for i in range(N_samples):
        base_sample_i = _get_single_clean_sample(sample_data, i, min_value=min_value)
        for j in range(i+1,N_samples):
            base_sample_j = _get_single_clean_sample(sample_data, j, min_value=min_value)
            # get rid of any points with 'nan' values in either of the two samples 
            # TODO if the data is growthrates, which may contain nan/inf/neginf, clean that up somehow!
            if verbose:
                print "Plotting sample %s (length %s) vs sample %s (length %s) (out of %s samples)"%(i,len(base_sample_i), 
                                                                                     j,len(base_sample_j), N_samples)
            sample_i, sample_j = remove_nan_inf_from_two_samples(base_sample_i, base_sample_j, min_val=min_val_to_plot)
            if not (sample_i and sample_j):
                print "Error: no datapoints left in samples %s and %s after removing nan/inf values!"%(i,j)
                continue
            if verbose:
                print "After removing nan/inf values and values below %s, %s datapoints left."%(min_val_to_plot, 
                                                                                                len(sample_i))
            mplt.subplot(N_samples, N_samples, N_samples*(i)+(j+1))
            if grey_out_min is not None:
                if len(set(colors)) > 1:
                    print "Warning: overwriting any colors set with -C/M/B/G/R options due to grey_out_min being set!"
                colors = [('k' if x>=grey_out_min and y>=grey_out_min else '0.6') for x,y in zip(sample_i,sample_j)]
            # note: i is on the y axis, j on the x axis, to make it fit the label below
            mplt.scatter(sample_j, sample_i, s=dotsize, c=colors, edgecolors='none')
            ### setting up the axis scale, ticks and tick labels - cosmetic stuff
            max_i = max(sample_i);  min_i = min(sample_i)
            max_j = max(sample_j);  min_j = min(sample_j)
            if plot_scale=='log': 
                mplt.loglog()
                #mplt.xscale(plot_scale); mplt.yscale(plot_scale)   # this is an alternative version to mplt.loglog()
                # I changed all values under 0.1 (including 0.1) to 0.1's, thus the limits and labels
                ymin,ymax = min_i/3, max_i*3
                xmin,xmax = min_j/3, max_j*3
                #text_pos_x, text_pos_y = min_value/5, min_i*3/5
                text_pos_x, text_pos_y = max_j*2, min_i*4/5
                # makes it so the ticks stay as they were but the only labels are <0.1, 1 and 100
                # TODO actually the .1 here should depend on min_value... 
                if minimize_tick_labels:
                    locs,labels = mplt.xticks();    mplt.xticks(locs,['','0','1','','100'])
                    locs,labels = mplt.yticks();    mplt.yticks(locs,['','0','1','','100'])
            else:
                #mplt.ylim(-max_i/5, max_i*1.1)
                #mplt.xlim(-max_j/5, max_j*1.1)
                range_i = max_i - min_i
                range_j = max_j - min_j
                ymin,ymax = min_i-range_i*0.1, max_i+range_i*0.1
                xmin,xmax = min_j-range_j*0.1, max_j+range_j*0.1
                #text_pos_x, text_pos_y = -max_j/6, -max_i/12
                text_pos_x, text_pos_y = max_j, min_i
                # makes it so the ticks stay as they were but only the first (0) and third are labeled, for simplicity
                if minimize_tick_labels:
                    locs,labels = mplt.xticks();    mplt.xticks(locs, ['',str(int(locs[1])),'',str(int(locs[3]))])
                    locs,labels = mplt.yticks();    mplt.yticks(locs, ['',str(int(locs[1])),'',str(int(locs[3]))])
            if same_scale_both_plots:
                xymin = min(xmin,ymin)
                xymax = max(xmax,ymax)
                mplt.ylim(xymin,xymax)
                mplt.xlim(xymin,xymax)
            else:
                mplt.ylim(ymin,ymax)
                mplt.xlim(xmin,xmax)
            # axis labels on the leftmost plot of each row only; sample name and #reads, and #mutants for counts
            if i==j-1: 
                sample_label = _make_sample_label_text(sample_names[i], sample_i, datatype, 
                                                       min_value=(0 if plot_scale=='lin' else min_value), 
                                                       end_spaces=5, end_newlines=4)
                mplt.ylabel(sample_label, rotation=0, fontweight='bold', fontsize='x-large')
            # putting the title on the first subplot instead of the whole figure is a cheat, but I don't know how else.
            if need_title:
                mplt.title(fig_title, x=-1, horizontalalignment='left',verticalalignment='bottom', fontsize='x-large')  
                need_title = False
            # print Spearman or Pearson correlation coefficient on each graph (remove all nan/inf value positions first)
            if not print_correlation=='none':
                # TODO don't take the (min_value,min_value) mutants into account for correlations!  Or probably the mutants that are below the cutoff to print, either...
                if print_correlation=='spearman':   corr_coeff = spearmanr(sample_i,sample_j)[0]
                elif print_correlation=='pearson':  corr_coeff = corrcoef(sample_i,sample_j)[0][1]
                # the corr_coeff text color is based on value: <.25 black, .25-.50 blue, .50-.75 magenta, .75-1 red.
                corr_coeff_colors = ['k','b','m','r','r']
                color = corr_coeff_colors[int(abs(corr_coeff)*16)/ 4]
                mplt.text(text_pos_x, text_pos_y, 'corr %.2f'%corr_coeff, color=color, 
                          horizontalalignment='right',verticalalignment='top')
    # need a label for the last column, below the last plot
    # TODO if there are just two samples (i.e. one plot), change the top_newlines/end_spaces/end_newlines values to 0 in both this xlabel and the ylabel above to make things look better (may want to abstract all this into a subfunction)
    sample_label = _make_sample_label_text(sample_names[-1], sample_j, datatype, 
                                           min_value=(0 if plot_scale=='lin' else min_value), 
                                           end_spaces=0, top_newlines=2)
    mplt.xlabel(sample_label, rotation=0, fontweight='bold', fontsize='x-large', multialignment='right')
    return fig

# TODO another kind of correlation plot, single for all samples - make the x axis be mutants (ordered by various things - median abundance, abundance in a particular dataset, even location) and the y axis be abundance, with a color for each dataset, so each mutant would have N dots, one per dataset, and we could see how well the datasets cluster/deviate/etc.

def savefig(fig,figname,figtype='png'):
    """ Save fig (should be a matplotlib figure object) as file named figname.figtype. """
    mplt.figure(fig.number) # activate fig to save and close it
    # TODO dpi should probably be a command-line option!  Careful - higher size/dpi takes a LOT of memory!
    # 
    mplt.savefig(figname+'.'+figtype,bbox_inches='tight',pad_inches=0.5, dpi=300)
    mplt.close()


def read_joint_mutant_file(infile, which_read_subtype):
    """ Read file generated by mutant_join_datasets.py; return a readcount/growthrate table, sample list, and mutant 
    position/gene data table."""
    all_data = []
    for line in open(infile):
        # ignore comment lines, EXCEPT the header line
        if line.startswith('#') and not line.startswith('# chromosome'):    
            continue
        # ignore special lines in my test data output reference format
        if line.startswith('<REGEX>#') or line.startswith('<IGNORE>'):
            continue
        all_data.append(line.strip('\n').split('\t'))
    # count the number of samples based on header line
    valid_sample_header_prefixes = ['reads_in_', 'readcount_in_', 'total_in_', 'perfect_in_', 
                                    'growthrate_in_', 'TD_in_', 'TD_change_in_']
    N_sample_columns = len([x for x in all_data[0] if any([x.startswith(s) for s in valid_sample_header_prefixes])])
    if not N_sample_columns:
        sys.exit("No samples found in the data! Wrong input file format? Header:\n\"%s\""%all_data[0])
    # mutant_data is the general data (position/gene); sample_data is the sample readcounts/growthrates after that 
    #  (variable number of columns); annotation_data is everything after that.
    mutant_data = [[all_data[x][i] for i in range(8)] for x in range(len(all_data))]
    sample_data = [[all_data[x][i] for i in range(8,8+N_sample_columns)] for x in range(len(all_data))]
    annotation_data = [[all_data[x][i] for i in range(8+N_sample_columns,len(all_data[0]))] for x in range(len(all_data))]
    # separate header line from data in each
    mutant_data_headers, mutant_data         = mutant_data[0], mutant_data[1:]
    full_sample_headers, full_sample_data    = sample_data[0],  sample_data[1:]
    annotation_data_headers, annotation_data = annotation_data[0], annotation_data[1:]
    # now take just the perfect or imperfect readcounts if desired
    if which_read_subtype=='perfect':
        columns_to_keep = [n for n,val in enumerate(full_sample_headers) if 'perfect' in val]
        sample_headers = [x for n,x in enumerate(full_sample_headers) if n in columns_to_keep]
        sample_data = [[float(data_row[x]) for x in columns_to_keep] for data_row in full_sample_data]
    elif which_read_subtype=='imperfect':
        total_columns = [n for n,val in enumerate(full_sample_headers) if 'reads' in val]
        perfect_columns = [n for n,val in enumerate(full_sample_headers) if 'perfect' in val]
        assert len(total_columns)==len(perfect_columns)
        sample_headers = [x.replace('total','imperfect') for n,x in enumerate(full_sample_headers) if n in total_columns]
        sample_data = [[(float(data_row[t])-float(data_row[p])) for t,p in zip(total_columns,perfect_columns)] 
                       for data_row in full_sample_data]
    else:
        columns_to_keep = [n for n,val in enumerate(full_sample_headers) if not 'perfect' in val]
        sample_headers = [x for n,x in enumerate(full_sample_headers) if n in columns_to_keep]
        sample_data = [[x for n,x in enumerate(data_row) if n in columns_to_keep] for data_row in full_sample_data]
        sample_data = [[float(data_row[x]) for x in columns_to_keep] for data_row in full_sample_data]
    # grab just the sample names instead of full headers
    sample_names = [x.split('_in_',1)[-1] for x in sample_headers]
    #print "all_data: %s x %s"%(len(all_data), len(all_data[0]))
    #print "mutant_data: %s x %s"%(len(mutant_data), len(mutant_data[0]))
    #print "sample_data: %s x %s"%(len(sample_data), len(sample_data[0]))
    #print "annotation_data: %s x %s"%(len(annotation_data), len(annotation_data[0]))
    return sample_data, sample_names, mutant_data


def remove_mutants(all_data, file_with_mutants_to_remove, min_readcount_for_removal):
    """ ______________
    Doesn't modify the input."""
    # TODO write docstring, dammit!
    # TODO write somewhere on the plot that mutants were removed?
    sample_data, sample_names, mutant_data = all_data
    # get a list of (chrom,strand,min_pos) tuples for the mutants to remove from the dataset
    dataset_to_remove = Insertional_mutant_pool_dataset()
    dataset_to_remove.read_data_from_file(file_with_mutants_to_remove)
    mutants_to_remove = set([(mutant.position.chromosome, mutant.position.strand, str(mutant.position.min_position)) 
                             for mutant in dataset_to_remove])
    # go over all mutants in full dataset, only copy ones not in mutants_to_remove to the new datasets
    new_sample_data, new_mutant_data = [], []
    for i in range(len(sample_data)):
        chrom,strand,min_pos = mutant_data[i][:3]
        if not (chrom,strand,min_pos) in mutants_to_remove:
            new_mutant_data.append(mutant_data[i])
            new_sample_data.append(sample_data[i])
    # remove samples that now have no reads
    new_sample_names = sample_names
    # making this a reverse range because otherwise removing earlier columns would mess up later ones
    for i in range(len(new_sample_names)-1,-1,-1):
        total_reads = sum([float(data_line[i]) for data_line in new_sample_data])
        if total_reads==0:
            print "Removing sample %s from data due to 0 total reads!"%new_sample_names[i]
            del new_sample_names[i]
            for x in range(len(new_sample_data)):
                del new_sample_data[x][i]
    return new_sample_data, new_sample_names, new_mutant_data


def get_mutant_colors(color_cassette_mutants=False, color_mutants_from_files={}, read_count_cutoff=1, grey_color='0.6'):
    """ Return a (chrom,strand,min_pos):color dictionary (default black) for plotting mutants. 
    If color_cassette_mutants is True, mutants with insertion_cassette as chromosome will be grey (0.6).
    color_mutants_from_files is a color:file dict - any mutant with at least read_count_cutoff reads 
     in one of the files will be assigned that color (what happens when a mutant is in multiple files is undefined). 
     (if color_cassette_mutants is True, cassette mutants will be grey regardless of presence in files).
    """
    # TODO write somewhere on the plot what the colors mean!
    if color_cassette_mutants:
        default_color_function = lambda pos_data: grey_color if is_cassette_chromosome(pos_data[0]) else 'black'
        if verbose:
            print "insertion_cassette mutants will be colored grey (%s)"%(grey_color)
    else:
        default_color_function = lambda pos_data: 'black'
    color_dict = general_utilities.keybased_defaultdict(default_color_function)
    # TODO what's the color precedence order here?  Currently undefined.
    for (color,mutant_file) in color_mutants_from_files.items():
        if verbose:
            print "mutants with %s+ reads in file %s will be colored %s"%(read_count_cutoff, mutant_file, color)
        dataset = Insertional_mutant_pool_dataset()
        dataset.read_data_from_file(mutant_file)
        for mutant in dataset:
            chrom, strand, min_pos = mutant.position.chromosome, mutant.position.strand, mutant.position.min_position
            # if we're coloring cassette mutants grey, exclude them from this
            if color_cassette_mutants and chrom=='insertion_cassette':  
                continue
            if mutant.total_read_count >= read_count_cutoff: 
                color_dict[(chrom,strand,min_pos)] = color
                color_dict[(chrom,strand,str(min_pos))] = color
    if verbose and color_mutants_from_files:
        print "%s colored mutants from %s files"%(len(color_dict),len(color_mutants_from_files))
    assert color_dict[('chromosome_1','+',100)] == 'black'
    assert color_dict[('insertion_cassette','+',100)] == grey_color if color_cassette_mutants else 'black'
    return color_dict


def main_functionality(infile, options=None):
    """ Main functionality - see module docstring for more info."""
    # if no options given, use defaults given by the parser
    if options is None:
        (options, _) = define_option_parser()("")

    # instead of passing the verbose option around to each function, just make it global (TODO is that really a good idea?)
    global verbose
    verbose = options.verbose

    # TODO implement -g option!  Or remove it if this is already done another way.
    if options.growthrate_infiles:      sys.exit("-g/--growthrate_infiles option NOT IMPLEMENTED! Exiting.")
    # parse the options for coloring mutants from specific files into a single dictionary
    color_mutants_from_files = {}
    if options.color_mutants_magenta:   color_mutants_from_files['magenta'] = options.color_mutants_magenta
    if options.color_mutants_blue:      color_mutants_from_files['blue'] = options.color_mutants_blue
    if options.color_mutants_green:     color_mutants_from_files['green'] = options.color_mutants_green
    if options.color_mutants_red:       color_mutants_from_files['red'] = options.color_mutants_red

    # read input data (all_data will be an (sample_data, sample_headers, mutant_data) tuple)
    all_data = read_joint_mutant_file(infile, options.which_reads)
    if options.remove_mutants_from_file:
        all_data = remove_mutants(all_data, options.remove_mutants_from_file, options.remove_mutants_readcount_min)
    # TODO it would be nice to be able to give more detailed names/descriptions for the datasets... How?

    # TODO right now this just takes a single joint-mutant file (generated by mutant_join_datasets.py as input - make it also take any number of single-mutant files (generated by mutant_count_alignments.py)
    # deal with the output filename/format arguments
    if options.outfile_prefix:  options.outfile_prefix += '_'
    options.outfile_prefix = os.path.join(options.outfile_directory,options.outfile_prefix)
    if options.graph_scale=="lin":     graph_scales = [('linear','_lin')]
    elif options.graph_scale=="log":   graph_scales = [('log','_log')]
    elif options.graph_scale=="both":  graph_scales = [('linear','_lin'),('log','_log')]
    # make a mutant-to-color dictionary from color_cassette_mutants and color_mutants_from_files
    mutant_to_color = get_mutant_colors(options.color_cassette_mutants, color_mutants_from_files)
    # turn off interactive mode to prevent the images showing on the screen and eating up memory
    mplt.ioff() 
    # for all files, make correlation and/or distribution plots (linear and/or logscale)
    # TODO doing all the plotting twice for linear and log scale is slow and stupid...
    # TODO this is ridiculously slow sometimes!  Doing it with several files on the whole-genome library takes up a lot of memory, apparently.  Might be better to rewrite it to only read the files it's currently using, one or two at a time?  It'll be slower but use less memory.  (make it an option, maybe?)
    for (figscale,suffix) in graph_scales:
        if not options.graph_type=='dist':
            figname = options.outfile_prefix + ('growthrate' if options.growthrate_infiles 
                                                else 'readcount')
            fig = plot_all_correlations(all_data, figscale, figname, 
                                        options.print_correlations, mutants_to_color=mutant_to_color)
            savefig(fig, figname+suffix, options.graph_format)
        if not options.graph_type=='corr':
            figname = options.outfile_prefix + ('growthrate_distributions' if options.growthrate_infiles 
                                                else 'readcount_distributions')
            value_type = 'growth rate' if options.growthrate_infiles else 'read count'
            fig = plot_all_distributions_sorted(all_data, figscale, figname, value_type)
            # TODO do I want mutant_to_color to apply for the distribution plot too?
            savefig(fig, figname+suffix, options.graph_format)
    # the return value is only for interactive use
    return all_data

# TODO for correlation plots, also have a variant where instead of passing one infile list A,B,C,D and getting an all-by-all plot, you pass two lists A,B C,D and just get a 2x2 (or 2x3, etc) correlation plot of first list vs second list?  Sometimes that's all you want, especially if you want the correlation of just one sample against several.


def define_option_parser():
    """ Returns the optparse OptionParser object - run "(options, infiles) = parser.parse_args()" to use. """
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-g', '--growthrate_infiles', action='store_true', default=False, 
                      help='interpret the numbers in infiles as growth rates (default: as deepseq read counts)')
    parser.add_option('-o', '--outfile_prefix', default='', help='the prefix to add to all output file names')
    parser.add_option('-d', '--outfile_directory', default='.', help='directory to create the output files in (default .)')
    parser.add_option('-t', '--graph_type', choices=['dist','corr','both'], default='corr', metavar='[dist|corr|both]')
    parser.add_option('-s', '--graph_scale', choices=['lin','log','both'], default='log', metavar='[lin|log|both]')
    parser.add_option('-f', '--graph_format', default='png', metavar='FORMAT', 
                      help='Which formats work depends on matplotlib backend - try a few and see. Default %default.')
    parser.add_option('-c', '--print_correlations', choices=['none','pearson','spearman'], default='spearman', 
                      metavar='none|spearman|pearson', 
                      help="Calculate and print the data correlation for each correlation graph (default %default).")
    parser.set_defaults(which_reads='all')
    parser.add_option('-p', '--perfect_reads_only',   action='store_const', dest='which_reads', const='perfect', 
                      help="Take only perfectly aligned read counts (default %default)")
    parser.add_option('-i', '--imperfect_reads_only', action='store_const', dest='which_reads', const='imperfect',
                      help="Take only imperfectly aligned read counts (default %default)")
    parser.add_option('-a', '--all_reads',            action='store_const', dest='which_reads', const='all',
                      help="Take total read counts (turns -p/-i off).")
    parser.add_option('-X', '--remove_mutants_from_file', metavar='FILE',
                      help='Remove all mutants present in FILE from the datasets (see -Y for read count cutoff).')
    parser.add_option('-Y', '--remove_mutants_readcount_min', type='int', default=1, metavar='N',
                      help='When applying -X, only remove mutants with at least N reads in FILE (default %default).')
    # TODO also add options to remove cassette mutants?  Or not, that can be done at the mutant_count_alignments.py stage

    parser.add_option('-v', '--verbose', action='store_true', default=False)

    # color options - TODO these should have a separate group, with stuff explained just once on top!  And specify the order of precedence (i.e. if a mutant is in the -M and -R files, will it be blue or red?
    # TODO all these should actually take multiple files!
    # TODO let the user specify the precedence order of the colors?
    parser.add_option('-C', '--color_cassette_mutants', action='store_true', default=False, 
                      help='color cassette mutants on the plots grey (default %default).')
    parser.add_option('-M', '--color_mutants_magenta', metavar='FILE',
                      help='color mutants present in FILE magenta on the plots (default %default).')
    parser.add_option('-B', '--color_mutants_blue', metavar='FILE',
                      help='color mutants present in FILE blue on the plots (default %default).')
    parser.add_option('-G', '--color_mutants_green', metavar='FILE',
                      help='color mutants present in FILE green on the plots (default %default).')
    parser.add_option('-R', '--color_mutants_red', metavar='FILE',
                      help='color mutants present in FILE red on the plots (default %default).')
    # TODO how should the coloring work for the readcount distribution plots, where each sample gets a color? Either grey out all cassette dots, or make them paler or something... 
    # MAYBE-TODO There should also be a -G option for other mutants to be greyed-out instead of or in addition to the cassette; other colors should be ignored.
    return parser


if __name__ == '__main__':
    parser = define_option_parser()
    (options, infiles) = parser.parse_args()
    try:
        [infile] = infiles
    except ValueError:
        parser.print_help()
        parser.error("Error: exactly one input file required!")
    # main functionality done by "main_functionality" function
    main_functionality(infile, options)

