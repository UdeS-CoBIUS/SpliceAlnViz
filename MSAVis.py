#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""


moduleauthor::  Safa Jammali
Université de Sherbrooke, Canada
for questions email us at Safa.Jammali@USherbrooke.ca

2018

"""

# import external bilio
import dash 
import flask
import io
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

#from plotStrucuture import *
import plotly.graph_objs as go
from plotly import tools
from plotly.grid_objs import Column
from Bio import SeqIO
import skbio
from skbio.alignment import global_pairwise_align_nucleotide
from skbio.sequence import DNA
import numpy as np
# for creating a graph object
import plotly.figure_factory as ff
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from skbio import DNA, TabularMSA
import os
from Bio import AlignIO
import time
import pandas as pd
from pandas import ExcelWriter
import urllib


import xlwt 
from xlwt import Workbook 


EVALUE = 1.0/10000000


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

#app.config['suppress_callback_exceptions']=True

##########real input files#################################################

def read_MSA_file(input_MSA_file):
	"""
	    This function allows to read a multiple alignment from a fasta file

	    Parameters
	    ----------

	    input_MSA_file: string
			MSA file path

	    Returns
	    -------
	    dict_msa

	    """
	id_seq=[]
	sequence=[]
	dict_msa={}
	for record in SeqIO.parse(input_MSA_file, "fasta"):
		id_seq.append(record.id)
		sequence.append(str(record.seq))

		dict_msa[record.id]=record.seq
	msa = skbio.alignment.TabularMSA.read(input_MSA_file, constructor=DNA)

	#return id_seq, sequence, msa, dict_msa
	return dict_msa

def read_structure_file(input_macro_file):

    """
	    This function allows to read sequences exon positions from txt file

	    Parameters
	    ----------

	    input_macro_file: string
			structure file path

	    Returns
	    -------
	    exon_dict:

	    """
    
    exon_dict = dict()
    j=0
    f = open(input_macro_file.encode('utf-8'))
    for line in f:
        line = line.strip()
       
        
        if len(line) > 3:
            
            
            id_seq = line.split(":")[0]
            
            exon = line.split(":")[1]
            
            if id_seq  not in exon_dict:
                exon_dict[id_seq ] = []
            exon_dict[id_seq ].append(exon.split("-"))
   
    return exon_dict

def read_gene_cds(input_gene_cds_file):
	"""
	    This function allows to read the belonging cds to gene file

	    Parameters
	    ----------

	    input_gene_cds_file: string
			belonging cds to gene file path

	    Returns
	    -------
	    dict_seq_to_gene_belonging: dictionnary
			dict of gene (as key) and cds as value

	    """
	dict_seq_to_gene_belonging={}
	f = open(input_gene_cds_file.encode('utf-8'))
	for line in f:
		line = line.strip()
		[cds_id, gene_id]=line.split("\n")[0].split(" ")

		if len(dict_seq_to_gene_belonging)==0:
			dict_seq_to_gene_belonging[gene_id] = []
			dict_seq_to_gene_belonging[gene_id].append(cds_id)
		elif gene_id not in dict_seq_to_gene_belonging.keys():
			dict_seq_to_gene_belonging[gene_id] = []
			dict_seq_to_gene_belonging[gene_id].append(cds_id)
		else :
			dict_seq_to_gene_belonging[gene_id].append(cds_id)



	return dict_seq_to_gene_belonging

def read_fasta_file(file__sequence):
	dict_seq={}

	id_seq=[]
	sequence=[]
	
	for record in SeqIO.parse(file__sequence, "fasta"):
		#id_seq.append(record.id)
		#sequence.append(str(record.seq))
		dict_seq[record.id]=record.seq
	return dict_seq
##################################################################2 nd function############################################################
def compute_predicted_exon(dict_msa, dict_seq_to_gene_belonging, exon_position ):
	predicted_exon={}
	
	listGene=dict_seq_to_gene_belonging.keys()
	listCDS=dict_seq_to_gene_belonging.values()
	dict_aln_to_plot= convert_dict_position_to_MSA_pos(dict_msa, exon_position, listGene,listCDS )
	for geneid, listCDS in dict_seq_to_gene_belonging.items():
		
		predicted_exon[geneid]=[]
		
		list_exon_gene=dict_aln_to_plot[geneid]
		
		list_exon_CDS=[]
		for idCDS in listCDS:
			predicted_exon[idCDS]=[]
			for i in dict_aln_to_plot[idCDS]:
				list_exon_CDS.append(i)
			
		#list_exon_CDS=list(set(list_exon_CDS))
		for i in list_exon_gene:
			if i not in list_exon_CDS:
				for j in range (i[0], i[1]+1):
					predicted_exon[geneid].append(j)
	
	return predicted_exon
def plotMSA(height_value, width_value,input_MSA_file,input_macro_file, input_gene_cds_file,  select_features, chosen_sequences, search_motif ):


        
	#read input files
	exon_position= read_structure_file(input_macro_file.encode( 'utf-8'))
	dict_msa= read_MSA_file(input_MSA_file.encode('utf-8'))
	dict_seq_to_gene_belonging = read_gene_cds(input_gene_cds_file)
	
	#prepare positions
	start_codon, stop_codon = search_position(dict_msa)
	
        listGene=dict_seq_to_gene_belonging.keys()
	listCDS=dict_seq_to_gene_belonging.values()
    	flat_list_listCDS=[]
	for sublist in listCDS:
	    for item in sublist:
               	flat_list_listCDS.append(item)
	list_CDS=list(set(flat_list_listCDS))
	
	splice_sites_pos= extrat_splice_sites_position(exon_position, listGene,chosen_sequences )
	
	splice_sites= convert_postion_to_MSA_pos(dict_msa, splice_sites_pos)
	dict_intron={}
	
	
	exon_junction_pos= extrat_exon_junction_position(exon_position, list_CDS,chosen_sequences )
	exon_junction= convert_postion_to_MSA_pos(dict_msa, exon_junction_pos)

	predicted_exon=compute_predicted_exon(dict_msa, dict_seq_to_gene_belonging, exon_position )

	dict_choosen_sequence, base_text, base_values, list_order_key= msa_defaultplot(chosen_sequences, dict_msa)
	
	#select_features==['predicted exon']:

	#elif select_features == [ 'frame1']:
	#	position_to_display=frame1_position
	#elif select_features == [ 'frame2']:
	#	position_to_display=frame2_position
	#elif select_features == [ 'frame3']:
	#	position_to_display=frame3_position
	
	#elif select_features == [ 'frameshift']:
	#	position_to_display=frameshift_position
	#elif select_features == [ 'intron']:
	#	position_to_display=intron_position
	#elif select_features == [ 'splice sites']:
	#	position_to_display=splicesites_position
	if  'start codon' in select_features:
		
		base_values = hmap_toplot_modif(start_codon,dict_choosen_sequence, base_values, list_order_key, 'codon start')
	
	if 'stop codon' in select_features:
		
		base_values = hmap_toplot_modif(stop_codon,dict_choosen_sequence, base_values, list_order_key, 'codon stop')
	if 'splice sites' in select_features:
		
		base_values = hmap_toplot_modif(splice_sites,dict_choosen_sequence, base_values, list_order_key, 'splice sites')
	if 'exon junction' in select_features:
		
		base_values = hmap_toplot_modif(exon_junction,dict_choosen_sequence, base_values, list_order_key, 'exon junction')
	
	if 'predicted exon' in select_features:
		
		base_values = hmap_toplot_modif(predicted_exon,dict_choosen_sequence, base_values, list_order_key, 'predicted exon')
	
		
	id_sequence_to_viz = list_order_key
	
	
	if search_motif:
		print 'enter',search_motif
	#fig_ = ff.create_annotated_heatmap(base_values, annotation_text=base_text, colorscale=colorscale)
	flat_list=[]
       
	for sublist in base_text:
		for item in sublist:
			flat_list.append(item)
        annotations = go.Annotations()
	
        msa_consensus=[]
	list_sequence_to_display=[]
	for bloc in dict_choosen_sequence.values():
                
		msa_consensus.append(DNA(str(bloc[1])))
	
        	list_sequence_to_display.append(str(bloc[1]))

	consensus= find_consensus( TabularMSA(msa_consensus))
				
	conservation= compute_consevation(list_sequence_to_display)
	
	
	traceA, traceC, traceG, traceT, yA, yC, yG, yT=plot_consensus(consensus, conservation)
	
	fig = tools.make_subplots(rows=5, cols=1)

	
	consensus_=[str(consensus[i]) for i in range(len(consensus))]
	

	cordx= [i for i in range(len(consensus_))]
	cordy= [i for i in range(len(id_sequence_to_viz))]
	for n, row in enumerate(base_text):
	    
	    for m, val in enumerate(row):
		
		annotations.append(go.Annotation(text=str(base_text[n][m]), x=cordx[m], y=cordy[n],
                                 xref='x1', yref='y1', font=dict(family='Courier New',
				size=18,
				color='black',
				), showarrow=False))
	
	fig_ =  go.Heatmap(
		x=cordx, y=cordy, z=base_values.tolist(),
		showscale= False,
				
		colorscale=colorscale,
		
	    )
	
	fig.append_trace(fig_, 1, 1)
	fig.append_trace(traceA, 2, 1)
	fig.append_trace(traceC, 3, 1)
	fig.append_trace(traceG, 4, 1)
	fig.append_trace(traceT, 5, 1)
	
	height_value_= 100+( int(height_value)*int(len(id_sequence_to_viz)))
	width_value_= 100+ (int(width_value)*int(len(consensus_)))
	#for i in range(len(fig.layout.annotations)):
	 #   fig.layout.annotations[i].font.size = 18
	fig['layout'].update(
		annotations=annotations,

		xaxis=dict(
			#autorange=True,
			#showgrid=False,
			#zeroline=False,
			#showline=False,
			ticktext= consensus_,
			tickvals=[i for i in range(len(consensus_))],
			
			tickfont = {'family': 'Courier New','color':'red', 'size': 24}
			),

				
		
		yaxis=dict(autorange='reversed',
		   ticks='',
		   ticksuffix=' ',
		   ticktext=id_sequence_to_viz,
		   tickvals=list(np.arange(0, len(base_text[0]))),
		   showticklabels=True,
		   tickfont=dict(family='Courier New',
				size=18,
				color='#22293B',
				),
			
			domain=[0.5, 1],
		
			
		    ),
		#width=1000,
		#height=450,
		#
		#layout['annotations'] = annotations	
	
           
           	
		xaxis2=dict(
			#autorange=True,
			#showgrid=False,
			#zeroline=False,
			showline=False,
			#ticktext=consensus,
			tickvals=[i for i in range(len(yA))],
			ticks='',
			
			showticklabels=False),
		
		xaxis3=dict(
			#autorange=True,
			#showgrid=False,
			#zeroline=False,
			showline=False,
			#ticktext=consensus,
			tickvals=[i for i in range(len(yA))],
			ticks='',
			showticklabels=False),
		
		xaxis4=dict(
			#autorange=True,
			#showgrid=False,
			#zeroline=False,
			showline=False,
			#ticktext=consensus,
			tickvals=[i for i in range(len(yA))],
			ticks='',
			showticklabels=False),
		xaxis5=dict(
			#autorange=True,
			#showgrid=False,
			#zeroline=False,
			showline=False,
			ticktext=[i for i in range(len(yA))],
			tickvals=[i for i in range(len(yA))],
			tickfont=dict(family='Courier New',
				size=18,
				color='#22293B',
				),
			#showticklabels=False
			),
		yaxis2=dict(title='A',
			titlefont=dict(
			    size=20),
			domain=[0.3, 0.4],
			color='orange',
			showticklabels=False
		    ),
		yaxis3=dict(title='C',
			titlefont=dict(
			    size=20),
			domain=[0.2, 0.3],
			color='green',
			showticklabels=False
		    ),
		yaxis4=dict(title='G',
			titlefont=dict(
			    size=20),
			domain=[0.1, 0.2],
			color='red',
			showticklabels=False
		    ),
		yaxis5=dict(title='T',
			titlefont=dict(
			    size=20),
			domain=[0, 0.1], 
			color='purple',
			
			showticklabels=False
		    ),
		#rangeslider=dict(
		#		    visible = True
		#		),
		autosize = False,
		height= height_value_,
		width= width_value_,
		showlegend = False,
		)
		
	
	
				
	return fig
def search_position(dict_msa):
	startcodon_position={}
	

	stopcodon_position={}
	

	for key_id_CDS in dict_msa.keys():
		startcodon_position[key_id_CDS]=[]
		stopcodon_position[key_id_CDS]=[]
		sequence_value=str(dict_msa[key_id_CDS])
		
		for position_nt in range (0, len(sequence_value)-2):
			codon=sequence_value[position_nt]+sequence_value[position_nt+1]+sequence_value[position_nt+2]
                        
			if (codon == start):
				startcodon_position[key_id_CDS].append(position_nt)
				startcodon_position[key_id_CDS].append(position_nt+1)
				startcodon_position[key_id_CDS].append(position_nt+2)
			if (codon in  stop_codon):
				stopcodon_position[key_id_CDS].append(position_nt)
				stopcodon_position[key_id_CDS].append(position_nt+1)
				stopcodon_position[key_id_CDS].append(position_nt+2)

	
	return startcodon_position, stopcodon_position
#exon_position_in_msa={'s1':[[1,5],[10, 15]]}

def extrat_splice_sites_position(exon_position, listGene,chosen_sequences ):
	splice_sites={}

	
	
	for key in chosen_sequences:
                if key in listGene:
			list_one_gene= exon_position[key]

			
			splice_sites_list=[]
		        #length=int(list_one_gene[0][1]) - int(list_one_gene[0][0])
			length=int(list_one_gene[0][1])-int(list_one_gene[0][0])
			splice_sites_list.append(length)
			
			for exon in list_one_gene[1:-1]:
				len_exon= int(exon[1])-int(exon[0])
				
				length=length+len_exon
				splice_sites_list.append(length)
				
			splice_sites[key]=splice_sites_list
		else:
			splice_sites[key]=[]
	
	return splice_sites

def extrat_exon_junction_position(exon_position, listCDS,chosen_sequences):
	exon_junction={}
	
	for key in chosen_sequences:
                if key in listCDS:
			list_one_CDS= exon_position[key]
			exon_junction_list=[]
		        
			for exon in list_one_CDS[1:]:
			
				exon_junction_list.append(int(exon[0]))
			exon_junction[key]=exon_junction_list
		else:
			exon_junction[key]=[]
	
	return exon_junction

def convert_postion_to_MSA_pos(dict_msa, exon_pos):

	
	new_dict_pos={}
	for key in exon_pos:
		
		if len (exon_pos[key])>0:
			
			new_pos=convert_listpostion_to_msaalignment(str(dict_msa[key]),exon_pos[key])
			new_dict_pos[key]= new_pos
		else:
			new_dict_pos[key]= []
	
	return new_dict_pos
def convert_listpostion_to_msaalignment(ref_sequence, listposition):
	"""
	    This function allows 

	    Parameters
	    ----------

	   

	    Returns
	    -------
	    


	    """
	sequence= ref_sequence.replace("-","")
        
	
	list_position_aln=[]
	
	
	cmpt=0
        count_element=0
	for ind, nucleo in enumerate(ref_sequence):
		
		
		if nucleo!='-':
			
			if cmpt in listposition:
				list_position_aln.append(ind)
				
			
				count_element=count_element+1
			
			cmpt = cmpt+1
	
	
	return list_position_aln

def msa_defaultplot(chosen_sequences, dict_msa):

	
	
	base_text=[]
	dict_choosen_sequence=select_choosen_sequence( chosen_sequences, dict_msa)
	
	list_order_key=[]
	for key in dict_choosen_sequence.keys():
		list_order_key.append(key)
		base_text.append(list(str(dict_choosen_sequence[key][1])) )
		
	
	base_values = np.zeros((len(base_text), len(base_text[0])))
		
	return(dict_choosen_sequence,base_text, base_values, list_order_key)
def select_choosen_sequence( list_seq, dict_msa):
	"""
		    This function allows to read a belonging cds to gene file

		    Parameters
		    ----------



		    Returns
		    -------


		    """

	dict_choosen_sequence={}
	
	for key in dict_msa.keys():
		if key in list_seq:

			dict_choosen_sequence[key]=[list_seq.index(key), dict_msa[key]]

	return dict_choosen_sequence

def hmap_toplot_modif(position_to_display, dict_choosen_sequence, base_values, list_order_key, component):
	
	
	for i in range (0, len(dict_choosen_sequence.keys())):
                        
			current_sequence_id= dict_choosen_sequence.keys()[i]
			for j in position_to_display[current_sequence_id]:
					
					base_values[i][j] = base_dic[component]

					 
						

	
	return base_values


def find_consensus(msa):
	"""
		    This function computes the consensus of a givin msa

		    Parameters
		    ----------

		    msa:
				alignment


		    Returns
		    -------
		    consensus:


		    """
	consensus= msa.consensus()
	return consensus

def compute_consevation(list_sequence):
	conservation=[]

	col=[]
	for i in range (0, len(list_sequence[0])):
		
		
		for j in range (0, len(list_sequence)):
			
			col.append(list_sequence[j][i])
	
		column_identity= percent_id_per_column(col)
		conservation.append(column_identity)
		col=[]
	return conservation
def percent_id_per_column(column):
	percent_identity=0
	l = [i for i in column if i != '-']
	if len(l)==0:
		percent_identity= [0,0,0,0]
	else:
		percent_identity = [0,0,0,0]
		A= column.count('A')
		C= column.count('C')
		G= column.count('G')
		T= column.count('T')
		percent_identity = [100*A/len(l),100*C/len(l),100*G/len(l),100*T/len(l)]
	return percent_identity 



def plot_consensus(consensus, conservation):
	

	x=[i for i in range(0,len(consensus))]
	yA=[] 
	yC=[]
	yG=[]
	yT=[]
	
	for col in conservation:
		
		yA.append(col[0]) 
		yC.append(col[1])
		yG.append(col[2])
		
		yT.append(col[3])
	

	traceA = go.Bar(
	  x= x,
	  y= yA,
	  text=yA,
          textposition = 'inside',
	  showlegend= False,
	  textfont=dict(family='Courier New',
				size=18,
				color='black',
				)
	 )
	
	traceC = go.Bar(
	  x= x,
	  y= yC,
	  text=yC,
          textposition = 'inside',
	  showlegend= False,
	  textfont=dict(family='Courier New',
				size=18,
				color='black',
				)
	 )

	
	traceG = go.Bar(
	  x= x,
	  y= yG,
	  text=yG,
          textposition = 'inside',
	  showlegend= False,
	  textfont=dict(family='Courier New',
				size=18,
				color='black',
				)
	 )

	traceT = go.Bar(
	  x= x,
	  y= yT,
          text=yT,
          textposition = 'inside',
	  showlegend= False,
	  textfont=dict(family='Courier New',
				size=18,
				color='black',
				)
	  )

	data = [traceA, traceC, traceG, traceT ]
	layout = {
	  #'xaxis': {'title': 'X axis'},
	  #'yaxis': {'title': 'Y axis'},
	  'barmode': 'relative',
	 # 'title': 'Relative Barmode'
	}
	#data = go.Bar(
		
	#	x=[i for i in range(0,len(consensus))],
	#	y=conservation,
	#	showlegend= False,

    
    	#	text=consensus,
	#	textposition='auto'
		#showlegend=False,
		#text=consensus,
		#text=[i for i in range(len(consensus))]
		
	#)


	'''
	layout = go.Layout(
		xaxis=go.layout.XAxis(
			 ticktext=consensus,
			tickvals=[i for i in range(len(consensus))])
	)
	#fig = go.Figure(data=data, layout=layout)
	'''
	return  traceA, traceC, traceG, traceT, yA, yC, yG, yT
##################################################################Search motif ###########################################################

def search_motif_MSA(value, chosen_sequences, msa):
	
	list_of_file=[]
	for key  in msa:
		if key in chosen_sequences:
			
			MotifSeqFile = os.getcwd() + '/sequences/query/' + value + '.fasta'
			MotifFile = open(MotifSeqFile, "w")
			MotifFile.write(">" + value + "\n")
			MotifFile.write(str(value))
			MotifFile.close()
			
			SeqFile = os.getcwd() + '/sequences/subject/' + key + '.fasta'
			SequenceFile = open(SeqFile, "w")
			SequenceFile.write(">" + key + "\n")
			subject_sequence= str(msa[key]).replace("-","")
			SequenceFile.write(subject_sequence)
			
			SequenceFile.close()
			
			blastoutput= launch_tblastx(MotifSeqFile,SeqFile,  value, key)
			list_of_file.append(blastoutput)
	
	return list_of_file
def launch_tblastx( MotifSeqFile,SeqFile,motifid, sequenceid):
	"""
	This function 
	Parameters
	----------

	cdsfile:
	genefile:
	evalue:
	block1:
	block2:
	cdsid:
	geneid:

	Returns
	-------
	blastoutput:
	"""
	
	
	blastoutput = str(os.getcwd()+'/sequences/results/blast_results/Motif_'+motifid+'_vs_'+sequenceid+'.xlsx')
    	command = "tblastx -query " + MotifSeqFile + " -subject " + SeqFile +  " -evalue " + str(EVALUE) + " -outfmt "+ str(7) +  " -strand plus " +  " -out " + str(blastoutput)
		
	os.system(command)
	return blastoutput     
##################################################################2 nd function############################################################
def extractandplotStructure(input_macro_file, input_MSA_file, input_gene_cds_file, file_gene_sequence, file_CDS_sequence, choice_list, height_value, width_value):

	dict_msa= read_MSA_file(input_MSA_file.encode('utf-8'))	
	exon_position= read_structure_file(input_macro_file.encode( 'utf-8'))
        dict_seq_to_gene_belonging = read_gene_cds(input_gene_cds_file.encode( 'utf-8'))
	gene_sequence= read_fasta_file(file_gene_sequence)
	CDS_sequence= read_fasta_file(file_CDS_sequence)


        listGene=dict_seq_to_gene_belonging.keys()
	listCDS=dict_seq_to_gene_belonging.values()
      	dict_aln_to_plot={}
	#for key in chosen_sequence:
	#	new_dict_msa[key]= dict_msa[key] 
	#	new_exon_pos[key]= exon_pos[key]
	
	dict_aln_to_plot= convert_dict_position_to_MSA_pos(dict_msa, exon_position, listGene,listCDS )
	listkey= dict_msa.keys()
	
        lenSeq=len(str(dict_msa[listkey[0]]))
	
	
	
	data, list_to_plot = data_static(dict_aln_to_plot, listGene, choice_list)
	
	maxlengthX= lenSeq # to add it
	maxlengthY= 2*(len(listGene)+len(listCDS))
	
   	height_value_= 100+ int(height_value)*int(maxlengthX)
	width_value_= 100+ int(width_value)*int(maxlengthY)
	
	layout = {
		'xaxis': {
			'range': [0, maxlengthX],
			'showgrid': False,
			 'showline': False,
			'zeroline':False,
		'tickvals':list(np.arange(0, maxlengthX)),
		'tickfont':dict(
				    family='Courier New',
				    size=18,
				    color='black'
				),
			
		    },
		'yaxis': {
			'range': [0, maxlengthY],
			 
        		#'zeroline':False,
        
			 'ticks':'',
			   'ticksuffix':' ',
			'showgrid': False,
			 'showline':False,
			'showticklabels':False
			
		    },
    		'shapes':list_to_plot,
		 'autosize':False,
		'height': height_value_,
		'width': width_value_,
		 
		}

	
	fig_static = go.Figure(data=data, layout=layout)
	
	return fig_static





def convert_dict_position_to_MSA_pos(dict_msa, exon_position, listGene,listCDS):
	"""
	    This function allows to convert coor
	    Parameters
	    ----------


	    Returns
	    -------
	 


	    """
	

	flat_listCDS=[]
		
	for sublistCDS in listCDS:
	    for item in sublistCDS:
       
		flat_listCDS.append(item)
       	
	new_dict_position={}
	#list gene
	for key in dict_msa.keys():
			
		
		if key in listGene:
			
			
			if len(exon_position[key])>0:
				
				msa_sites_list=[]
				begin=[]
				end=[]
				if len(exon_position[key])>1:
				


					
					
					length=0
					
					for exon in exon_position[key]:
						len_exon= int(exon[1])-int(exon[0])
				
						
						msa_sites_list.append([length,length+len_exon])
						length=length+len_exon


					
					for exon in msa_sites_list:
						begin.append(int(exon[0]))
						end.append(int(exon[1])-1)
				else:
					extrremity1=int(exon_position[key][0])
					extrremity2= extrremity1+(int(exon_position[key][1])-int(exon_position[key][0]))
					
					msa_sites_list.append([extrremity1,extrremity2])

					begin.append(extrremity1)
					end.append(extrremity2-1)
					
				
				pos_gene_begin= convert_listpostion_to_msaalignment(str(dict_msa[key]), begin)
				pos_gene_end= convert_listpostion_to_msaalignment(str(dict_msa[key]), end)
				
				list_pos_gene=[]
				for i in range (0, len(pos_gene_begin)):
					 
					list_pos_gene.append([int(pos_gene_begin[i]),int(pos_gene_end[i] )])
				
				new_dict_position[key]=list_pos_gene	
			else:
				new_dict_position[key]=[]
				

		
			
		#list CDS
		
		if key in flat_listCDS:
			
			if len(exon_position[key])>0:
				begin=[]
				end=[]	
				if len(exon_position[key])>1:
				
					
					for exon in exon_position[key]:
						begin.append(int(exon[0]))
						end.append(int(exon[1])-1)
				else:
					begin.append(int(exon_position[key][0]))
					end.append(int(exon_position[key][1])-1)
					
				
				pos_cds_begin= convert_listpostion_to_msaalignment(str(dict_msa[key]), begin)
				pos_cds_end= convert_listpostion_to_msaalignment(str(dict_msa[key]), end)
				
				list_pos_cds=[]
				for i in range (0, len(pos_cds_begin)):
					 
					list_pos_cds.append([int(pos_cds_begin[i]),int(pos_cds_end[i] )])
				
				new_dict_position[key]=list_pos_cds	
			else:
				new_dict_position[key]=[]
				
	return new_dict_position

	
'''
def convert_to_all_alignment_position(dict_choosen_sequence, dict_choosen_position):

	
	dict_junction={}
	for key in dict_choosen_sequence.keys():
		sequence=str(dict_choosen_sequence[key][1])
		listposition=dict_choosen_position[key]
		
		newposition =[]
		newposition= convert_to_alignment_position(sequence, listposition)
		dict_junction[key]= newposition
		
	return dict_junction


def hmap_toplot(position, chosen_sequences, dict_msa, component):

	list_sequence_reference=[]
	list_indice=[]
	base_text=[]

	dict_choosen_sequence=select_choosen_sequence( chosen_sequences, dict_msa)
	dict_choosen_position = select_choosen_position(chosen_sequences, position)
	
	for key in dict_choosen_sequence.keys():

		base_text.append(list(str(dict_choosen_sequence[key][1])) )
		
	#prepare 
	base_values = np.zeros((len(base_text), len(base_text[0])))

	
	#fill exon heatmap
	list_order_key=[]
	for i in range (0, len(dict_choosen_sequence.keys())):
			current_sequence_id= dict_choosen_sequence.keys()[i]
			list_order_key.append(current_sequence_id)
			#for j in range(0, len(str(dict_choosen_sequence[current_sequence_id][1]))):
				
			for bloc in dict_choosen_position[current_sequence_id]:
					
					for  j in range (int(bloc[0]) ,  int(bloc[1])):
						#i=dict_choosen_sequence[current_sequence_id][0]

						base_values[i][j] = base_dic[component]

					#else:
					#	base_values[i][j] = 
						

	#exit(-1)
	return(dict_choosen_sequence,base_text, base_values, list_order_key)

def select_choosen_position( list_seq, position):
	"""
		    This function allows to read a belonging cds to gene file

		    Parameters
		    ----------



		    Returns
		    -------


		    """

	dict_choosen_position={}

	for key in position.keys():
		if key in list_seq:

			dict_choosen_position[key]=position[key]

	return dict_choosen_position
'''



def data_static(dict_aln_to_plot, listGene, choice_list):
	
        list_to_plot=[]
	data=[]
	cmpt=0
	
	for key, val in dict_aln_to_plot.items():
		
		tracex= 'trace'+str(cmpt)
		tracex = go.Scatter(
		    #x=,
		    y= [1+cmpt],
	  	    text=[key],
		    textfont=dict(family='Courier New',
				size=18,
				color='black',
				),
		    #orientation = "v",
			# ['top left', 'top center', 'top right', 'middle left',
           # 'middle center', 'middle right', 'bottom left', 'bottom, center', 'bottom right']

		    textposition = 'bottom right',
	  	    showlegend= False,
		    mode='markers+text',
		    
		)
		data.append(tracex)
		if choice_list == 'CDS':
			
			if key not in listGene:
				
				for bloc in val:
					list_to_plot.append(creerShape(bloc[0],bloc[1] ,1+cmpt , 1.5+cmpt ,'rect', 'rgba(128, 0, 128, 0.7)'))
				
		elif choice_list == 'Gene':
			
			if key in listGene:
				for ind, bloc in enumerate(val):
					
					list_to_plot.append(creerShape(bloc[0],bloc[1] ,1+cmpt , 1.5+cmpt ,'rect', 'rgba(128, 0, 128,0)'))
					if len(val)>1 and ind< len(val)-1:
						list_to_plot.append(creerShape(bloc[1],val[ind+1][0] ,1+cmpt , 1+cmpt ,'line', 'rgba(128, 0, 128, 0)'))
				
		else:

			
			if key not in listGene:
				for bloc in val:
					list_to_plot.append(creerShape(bloc[0],bloc[1] ,1+cmpt , 1.5+cmpt ,'rect', 'rgba(128, 0, 128, 0.7)'))
			if key in listGene:
				for ind, bloc in enumerate(val):
					
					list_to_plot.append(creerShape(bloc[0],bloc[1] ,1+cmpt , 1.5+cmpt ,'rect', 'rgba(128, 0, 128, 0)'))
					if len(val)>1 and ind< len(val)-1:
						list_to_plot.append(creerShape(bloc[1],val[ind+1][0] ,1+cmpt , 1+cmpt ,'line', 'rgba(128, 0, 128, 0)'))
			

		cmpt= cmpt+1


       
	#fig = go.Figure(data=data, layout=layout)		
	#return fig
	return data,  list_to_plot


	
def creerShape(x0, x1,y0, y1, type_, fillcolor):
	if type_=='line':
		return  {  'type': type_,
			'x0': x0,
		    'y0': y0,
			'x1': x1,
		    'y1': y1,
		 #'fillcolor': 'rgba(128, 0, 128, 0.7)',
		}

	else:
		return  {  'type': type_,
			'x0': x0,
		    'y0': y0,
			'x1': x1,
		    'y1': y1,
		 'fillcolor': fillcolor,
		}
	



##############Display example######################################################################################


input_MSA_file=os.getcwd()+'/datas/msa_base.fasta'

input_macro_file= os.getcwd() +'/datas/structure_base.txt'
input_gene_cds_file= os.getcwd() + '/datas/source2Target_base.txt'
file_gene_sequence=os.getcwd() + '/datas/gene_base.fasta'
file_CDS_sequence=  os.getcwd()  +'/datas/cds_base.fasta'
#consensus=['G','C','C','C','T','A','G','G','A','G','A','A','C','G','C','G']
#conservation=[100,100,100,100,60,100,100,100,100,100,100,100,100,50,100,100]
#figshape={}
colorscale = [[0, '#F4F0E4'],#99CCCC
        [0.1, '#e7298a'],
	[0.2, '#1b9e77'],
	[0.3, '#d95f02'],
	[0.4, '#7570b3'],
	[0.5, '#e7298a'],
	[0.6, '#BA817A'],
	[0.7, '#DB7C1C'],
	[0.8, '#BFE231'],
	[0.9, '#179835'],
	[1, '#981717']
	]
base_dic = {'predicted exon': 0.1, 'intron_color': 0.2,'splice sites':0.3, 'exon junction': 0.4, 'codon start': 0.5,'codon stop': 0.6, 'frame1_color':0.7, 'frame2_color':0.8,'frame3_color':0.9,'motif_color':1  }
 
#component = 'exon'
#frame_position={'s2':[[3,6], [9,13]]}
#exon_position_in_msa={'s1':[[1,5],[10, 15]]}
select_features=[]
#dropdown_choice_heatmap='Sequence'
dropdown_choice='All'
search_motif=''
start='ATG'
stop_codon=['TAA','TAG','TGA']

height_value=80
width_value=100
height_value2=30
width_value2=130
dict_msa= read_MSA_file(input_MSA_file.encode('utf-8'))
chosen_sequences= dict_msa.keys()
id_seq= dict_msa.keys()
fig= plotMSA(height_value, width_value,input_MSA_file,input_macro_file, input_gene_cds_file,  select_features, chosen_sequences, search_motif )

fig_static= extractandplotStructure(input_macro_file.encode( 'utf-8'),input_MSA_file.encode('utf-8') , input_gene_cds_file.encode('utf-8'), file_gene_sequence.encode('utf-8'),file_CDS_sequence.encode('utf-8'), dropdown_choice, height_value2, width_value2)



##################################################################################################################
  
app.layout = html.Div([
	html.Br(),
	html.H1(children='SpliceAlnViz: Multiple Spliced Alignment visualiszation tool',style={'textAlign': 'center', 'color': '#bc80bd'}), #"#bebada","#fdb462","#fb8072","#d9d9d9","#bc80bd","#b3de69","#8dd3c7","#80b1d3","#fccde5","#ffffb3"
	html.Hr(),   
	html.Div([
	
		
        dcc.Upload(
		id='file_gene_sequence',
		children=html.Div([
		    'Gene Sequence ',
		    html.A('Select a File')
			],  style={ 'fontSize': 14})
			),
	dcc.Upload(
		id='file_CDS_sequence',
		children=html.Div([
		    'CDS Sequence ',
		    html.A('Select a File')
			],  style={ 'fontSize': 14})
			),
        dcc.Upload(
		id='filegenecds',
		children=html.Div([
			'Gene to CDS belonging  ',
			html.A('Select a File')
		], style={'fontSize': 14})

		),
	dcc.Upload(
		id='file-msa',
		children=html.Div([
		    'Multiple Sequence Alignment ',
		    html.A('Select a File')
			],  style={ 'fontSize': 14})
		),
       
	dcc.Upload(
		id='filestructure',
		children=html.Div([
		    'CDS position ',
		    html.A('Select a File')
		], style={ 'fontSize': 14})
		
    		),

	dcc.RadioItems(
		id='dropdown',
		options=[{'label': 'OK', 'value': 'OK'}, {'label': 'NO', 'value': 'NO'}],
		value='NO'
   		 ),
    
	
	#html.Br(),
    html.Label('Select features', style={ 'fontSize': 14}),
    dcc.Dropdown( id='select_features',style={ 'fontSize': 14, 'width':200},
        options=[
            {'label': 'predicted exon', 'value': 'predicted exon'},
            {'label': 'splice sites', 'value': 'splice sites'},
           # {'label': 'frame1', 'value': 'frame1'},
	    #{'label': 'frame2', 'value': 'frame2'},
            #{'label': 'frame3', 'value': 'frame3'},
	    {'label': 'start codon', 'value': 'start codon'},
	    {'label': 'stop codon', 'value': 'stop codon'},
	    #{'label': 'intron', 'value': 'intron'},
	    {'label': 'exon junction', 'value': 'exon junction'}
        ],
        #value=['predicted exon'],
        multi=True
    ),
   
    #html.Br(),
    #html.Label('Sequence or Artificial', style={ 'fontSize': 14}),
    #dcc.RadioItems(
     #   id='dropdown_choice_heatmap',
      #  options=[{'label': i, 'value': i} for i in ['Sequence', 'Artificial']],
       # value='Sequence'
    #),
    #html.Br(),
     

    html.Label('Display sequence', style={ 'fontSize': 14}),
    dcc.Checklist(id='chosen_sequences',
        options=[
		 {'label': i, 'value': i} for i in id_seq
            
        ],
        values=id_seq
	
    )

    	

	],	style={'columnCount': 1}),
  	
	 	html.Br(),
		 html.Label('Search Motif', style={ 'fontSize': 14}),
		 #dcc.Input(id='search_motif',value='', type='text'),
		html.Div(dcc.Input(id='search_motif', type='text')),
			html.Button('Submit', id='button'), 
			#html.Div(id='output-container-button',
			#children='Enter a value and press submit')
				
		html.A(
		'Download Data',
		id='output-container-button',
		download="",
		href="",
		target="_blank"
	   		 ),
	    	
	   	#dcc.Slider(
		#	min=0,
		#	max=9,
		#	marks={i: 'Label {}'.format(i) for i in range(10)},
		#	value=5,
		#	) ,
		html.Br()   , 
		html.Label('MSA Figure scale', style={ 'fontSize': 14}),
		dcc.Input(id='height_value',
			    placeholder='Enter a height value',
			    type='int',
			    #value=height_value
			),
		dcc.Input(id='width_value',
			    placeholder='Enter a width value',
			    type='int',
			    #value=width_value
			),
		html.Br(),
		html.Div([	
			dcc.Graph(id='output-data-upload',figure=fig) ], style={'overflowX': 'scroll', 'width': '1080'}),
				#, style={'marginBottom': 50, 'marginTop': 25}),
				
		
		
		#html.Br(),
		html.Div([	
			dcc.Graph(id='consensus_msa_output',figure=fig_static) 
				 
			#, style={'marginBottom': 50, 'marginTop': 25}
			], style={'overflowX': 'scroll', 'width': '1080'}),
		
		html.Label('MSA structure scale', style={ 'fontSize': 14}),
		dcc.Input(id='height_value2',
			    placeholder='Enter a height value',
			    type='int',
			    #value=height_value
			),
		dcc.Input(id='width_value2',
			    placeholder='Enter a width value',
			    type='int',
			    #value=width_value
			),
		 html.Label('Select sequence category', style={ 'fontSize': 14}),
		    dcc.RadioItems(
			id='dropdown_choice',
			options=[{'label': i, 'value': i} for i in ['All', 'CDS', 'Gene']],
			value='All'
		    ),
		html.Hr(),
		html.Label('Contact: Safa.Jammali@Usherbrooke.ca', style={ 'fontSize': 12,}),

],   className="container")
		 

@app.callback(
    dash.dependencies.Output('chosen_sequences', 'options'),
    [dash.dependencies.Input('dropdown', 'value'), dash.dependencies.Input('file-msa', 'contents')],

state=[dash.dependencies.State('file-msa', 'filename')]  )

def set_sequence_options(dropdown,filemsa_contents, filemsa):
  
    if filemsa is not None and dropdown == 'OK':
	
	dict_msa= read_MSA_file(filemsa.encode('utf-8'))
      
	return [{'label': i, 'value': i} for i in dict_msa.keys()]
    return None

@app.callback(
    dash.dependencies.Output('chosen_sequences', 'values'),
    [dash.dependencies.Input('chosen_sequences', 'options')])

def set_sequence_value(available_options):
    
    list_all=[]
    for i in available_options:
	list_all.append(i['value'])
    
    return list_all





@app.callback(dash.dependencies.Output('output-data-upload', 'figure'),
			  [dash.dependencies.Input('height_value', 'value'),
				dash.dependencies.Input('width_value', 'value'),
			dash.dependencies.Input('search_motif', 'value'),
			dash.dependencies.Input('chosen_sequences', 'values'),
				
		dash.dependencies.Input('select_features', 'value'),
	       #dash.dependencies.Input('dropdown_choice_heatmap', 'value'),
               dash.dependencies.Input('file-msa', 'contents'),
              dash.dependencies.Input('filestructure', 'contents'),
		dash.dependencies.Input('filegenecds', 'contents')],

              state=[dash.dependencies.State('file-msa', 'filename'),
              		dash.dependencies.State('filestructure', 'filename'),
			dash.dependencies.State('filegenecds', 'filename')]
                 )
def update_output(height_value, width_value, search_motif, chosen_sequences,select_features, 
			filemsa_contents, filestructure_contents, filegenecds_contents, filemsa, filestructure, filegenecds):

	
	
	if filemsa is not None and  filestructure is not None and filegenecds is not None:
		chosen_sequences=[x.encode('utf-8') for x in chosen_sequences]
	
		fig = plotMSA(height_value, width_value, filemsa.encode('utf-8'),  filestructure.encode( 'utf-8'),
				  filegenecds.encode('utf-8'), select_features, chosen_sequences, search_motif)
		
		    
		return fig

	return None


@app.callback(
    dash.dependencies.Output('output-container-button', 'href'),
    [dash.dependencies.Input('button', 'n_clicks'), 
	
	dash.dependencies.Input('chosen_sequences', 'values'),
	dash.dependencies.Input('file-msa', 'contents') ],

    [dash.dependencies.State('search_motif', 'value'),
	 dash.dependencies.State('file-msa', 'filename')])
	

def update_output(n_clicks, chosen_sequences,    filemsa_contents,value,filemsa):
	
	if filemsa is not None and  value is not None:
		 
		dict_msa= read_MSA_file(filemsa.encode('utf-8'))
		chosen_sequences=[x.encode('utf-8') for x in chosen_sequences]
		
		output_search_motif= search_motif_MSA(value, chosen_sequences, dict_msa)
		df_list = []
		lines= []
		for i in range (0, len(output_search_motif)):
			f = open(output_search_motif[i].encode('utf-8'))
			
			for line in f:
				lines.append(line.strip())
	
		df = pd.DataFrame(lines)
		
		csv_string = df.to_csv(index=False, encoding='utf-8')
		csv_string = "data:text/csv;charset=utf-8," + urllib.quote(csv_string)
			
		return csv_string		
    #return 'The input value was "{}" and the button has been clicked {} times'.format(
     #   value,
      #  n_clicks    )


    




@app.callback(dash.dependencies.Output('consensus_msa_output', 'figure'),
			  [
			dash.dependencies.Input('filestructure', 'contents'),
			   dash.dependencies.Input('file-msa', 'contents'),
			   dash.dependencies.Input('filegenecds', 'contents'),
				 dash.dependencies.Input('file_gene_sequence', 'contents'),
				 dash.dependencies.Input('file_CDS_sequence', 'contents'),
				dash.dependencies.Input('dropdown_choice', 'value'),
				dash.dependencies.Input('height_value2', 'value'),
				dash.dependencies.Input('width_value2', 'value') ],
		  state=[ dash.dependencies.State('filestructure', 'filename'),
			   dash.dependencies.State('file-msa', 'filename'),
			   dash.dependencies.State('filegenecds', 'filename'),
				dash.dependencies.State('file_gene_sequence', 'filename'),
				dash.dependencies.State('file_CDS_sequence', 'filename')]
                 )


#file-gene-sequence

def update_output_static( filestructure_contents,   filemsa_contents, filegenecds_contents,file_gene_sequence_content, file_CDS_sequence_content,dropdown_choice, height_value2, width_value2, filestructure, filemsa,
							  filegenecds, file_gene_sequence, file_CDS_sequence  ):
	if   filestructure is not None and  filemsa is not None and   filegenecds is not None:
		
		fig_static=extractandplotStructure(filestructure.encode( 'utf-8'),
				 filemsa.encode('utf-8'),filegenecds.encode( 'utf-8'), file_gene_sequence.encode('utf-8'), file_CDS_sequence.encode('utf-8'), dropdown_choice, height_value2, width_value2)
		    
		return fig_static

	return None

if __name__ == '__main__':
	app.run_server(debug=True)
