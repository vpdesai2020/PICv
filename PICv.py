
def __init_plugin__(app):
    app.menuBar.addmenuitem('Plugin', 'command',
        label='PICv',
        command=lambda: mytkdialog(app.root))

api_url = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/"
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import time
import re
import requests
import json
import csv
import pandas as pd
import sys
import numpy as np


import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.metrics import silhouette_score
import itertools
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import os.path
import pymol
from pymol import cmd

try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True
try:
    import tkFileDialog, tkMessageBox
except ImportError:
    import tkinter.filedialog
    import tkinter.messagebox


def open_file():
        global load_protein
        load_protein = tkFileDialog.askopenfilename(title='Open PDB File', filetypes=(("pdb files", "*.pdb"),))
        Label2.configure(text="Uploaded file : "+load_protein)
       
        
def check():
    
    if not len(Entry1.get()) == 0:
        uniprot_id_val = re.compile(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})')
        mo = uniprot_id_val.search(Entry1.get())
        if mo is not None:
            if not len(Entry2.get()) == 0:
                function()  
            else:
                text_label_2.configure(text='''Please enter Clusters''')
        else:
            text_label_1.configure(text="Please check Uniport ID")   
    else:
        text_label_1.configure(text="Please enter Uniport ID")
        
    
def function(*args, **kw):
    global json_response,X
    Label2.configure(text='''Running,...''')
    
    def get_protein_info(acc_id):
        
        url = api_url + acc_id
        time.sleep(0.01)
        response = requests.get(url)
        
        if response.status_code == 200:

            text_label_1.configure(text=acc_id+"  Id found")
            return json.loads(response.content.decode('utf-8'))
            
        else:
            return None
    
    acc_id = Entry1.get()

    
    json_response = get_protein_info(acc_id)
    if(json_response == None):
        
        text_label_1.configure(text='Invalid uniprot ID')
    
    def convert_json_reponse_to_dict(json_response):

        resultDict = {}
        index = 0

        for key in json_response.keys():

            for dataObj in json_response[key]['data']:

                accessionId = dataObj['accession']

                for residueObj in dataObj['residues']:
                    resId = residueObj['startIndex']

                    for pdbObj in residueObj['allPDBEntries']:
                        pdbId = pdbObj

                        resultDict[index] = {}

                        resultDict[index]['accessionId'] = accessionId
                        resultDict[index]['resId'] = resId
                        resultDict[index]['pdbId'] = pdbId

                        index+=1
        
        return resultDict
        
    resultDict = convert_json_reponse_to_dict(json_response)
    
    df = pd.DataFrame.from_dict(resultDict, orient='index')
    cols = df.columns.tolist()
    cols = cols[0:1] + cols[-1:] + cols[1:2]
    df = df[cols]

    group = df.groupby(['accessionId','pdbId'])
    df2 = group.apply(lambda x: x['resId'].unique())
    df2 = df2.apply(pd.Series)

    df2.values.sort()

    min_res_id = int(np.nanmin(df2.values))
    max_res_id = int(np.nanmax(df2.values))
    
    def generate_occurrence_matrix(frame):
    
        residue_numbers = range(min_res_id, max_res_id + 1)  
        
        cols = frame.columns[2:]
        data = frame[cols].values
        combined = data.flatten()
        
        unique_elements = np.unique(combined[~pd.isnull(combined)])
        output = [1 if res in unique_elements else 0 for res in residue_numbers]
        labels = [res for res in residue_numbers]
        
        return pd.Series(output, index = labels)
    
    df2.reset_index(inplace=True)
    df = df2.groupby(['accessionId']).apply(generate_occurrence_matrix)

    df2 = df.apply(lambda x: list(x.values), axis=1)

    proteins = np.array(df.index.tolist())
    X = df.values
    
    # def plot_elbow_coefficient_graph(range_n_clusters):
    #     sse = {}
    #     for k in range_n_clusters:
    #         kmeans = KMeans(n_clusters=k).fit(X)
    #         sse[k] = kmeans.inertia_ # Inertia: Sum of distances of samples to their closest cluster center
    #     plt.figure()
    #     plt.plot(list(sse.keys()), list(sse.values()))
    #     plt.xlabel("Number of clusters")
    #     plt.ylabel("SSE")
    #     plt.draw()
    #     plt.show(block=False)

    # Function to compute the optimal number of clusters based on silhouette scores
    def compute_optimal_k(range_n_clusters):
        global silhouette_score,optimal_k
        silhouette_scores_by_cluster_dict = dict()

        for n_clusters in range_n_clusters:
            clusterer = KMeans(n_clusters=n_clusters, random_state=1)
            cluster_labels = clusterer.fit_predict(X)
            # centers = clusterer.cluster_centers_

            # The silhouette_score gives the average value for all the samples.
            # This gives a perspective into the density and separation of the formed clusters
            silhouette_avg = silhouette_score(X, cluster_labels)

            # Add the silhouette score computed to a dictionary that maps cluster numbers to silhouette scores
            silhouette_scores_by_cluster_dict[n_clusters] = silhouette_avg
            print("For n_clusters =", n_clusters,
                    "The average silhouette_score is :", silhouette_avg)
            # Label2.configure(text="For n_clusters = "+str(n_clusters)+
                    # "The average silhouette_score is : "+str(silhouette_avg))

        # Find the optimal value of k (number of clusters) based on silhouette scores computed
        # The optimal number of clusters is that which has the highest associated silhouette score
        optimal_k = max(silhouette_scores_by_cluster_dict, key=silhouette_scores_by_cluster_dict.get)

        # Return the optimal number of clusters
        return optimal_k
        # -----------------------------------------------------
    

    #K Means Clustering
    k=int(Entry2.get())
    # k = int(input('Enter the number of clusters for performing K Means CLustering: 8'))
    kmeans = KMeans(n_clusters=k, random_state=0)
    kmeans.fit(X)

    cluster_map = pd.DataFrame()
    cluster_map['Protein'] = proteins
    cluster_map['Cluster'] = kmeans.labels_
    cluster_map.set_index('Protein', inplace = True)

    cluster_dict = {i: proteins[np.where(kmeans.labels_ == i)[0]] for i in range(kmeans.n_clusters)}
    #Dictionary of the cluster IDs and corresponding member protein IDs
    cluster_index_dict = {i: np.where(kmeans.labels_ == i)[0] for i in range(kmeans.n_clusters)}

    #Using tSNE for visualisation of clusters
    tsne = TSNE(n_components=2, verbose=1, init='pca', n_iter=5000, method='exact')
    reduced_data = tsne.fit_transform(np.concatenate((kmeans.cluster_centers_, X), axis = 0))

    cluster_centers = reduced_data[:kmeans.n_clusters, :]

    cluster_indices = cluster_dict.keys()
    x_coords = reduced_data[:len(cluster_indices),0]
    y_coords = reduced_data[:len(cluster_indices),1]
    
    
    

    # -----------elbow graph------------

    
    
    def plot_clusters():
        acc_id = Entry1.get()
        fig, ax = plt.subplots()

        for index, x, y in zip(cluster_indices, x_coords, y_coords):
            ax.annotate(index, (x, y))
        center_colors = np.random.uniform(low = 0.0, high = 0.7, size = (len(cluster_indices), 3))
        #ax.scatter(x_coords, y_coords, c = center_colors)
        colours = ListedColormap(["deepsalmon","deepblue","hotpink","chartreuse", "cyan", "darksalmon","dash","deepolive","deeppurple","firebrick"])

        point_colors = np.array([center_colors[cluster_map.loc[p].tolist()[0]] for p in proteins])
        ax.scatter(reduced_data[len(cluster_indices):, 0], reduced_data[len(cluster_indices):, 1], c = point_colors, marker = 'o',alpha=0.7,cmap=colours)
        # ax.legend()
        #plt.xlim(x_axis_min,x_axis_max)
        #plt.ylim(y_axis_min,y_axis_max)

        frame = plt.gca()
        frame.axes.get_xaxis().set_visible(False)
        frame.axes.get_yaxis().set_visible(False)
        plt.tick_params(
            axis='both',        # changes apply to the x-axis
            which='both',       # both major and minor ticks are affected
            bottom = False,     # ticks along the bottom edge are off
            top = False,        # ticks along the top edge are off
            left = False,       # ticks along the left edge are off
            right = False)      # ticks along the right edge are off
        plt.rcParams["figure.figsize"] = (10,10)

        #To save the figure
        
        plt.title(acc_id)
        # plt.draw()
        plt.savefig("Cluster_Plot.png")
        

        cluster_members = [', '.join(value) for value in cluster_dict.values()]
    
        ymd=time.strftime("_%Y%m%d-%H%M%S")
        FileName=acc_id+ymd+'.txt'
        for key, value in cluster_dict.items():
            # print('Cluster ', key, ' has proteins ', value)
            clu_val=('CLUSTER :  ', key, ' Has proteins --> ', value)
            d = ' '.join([str(elem) for elem in clu_val])
            with open(FileName, "a") as txt_file:
                for line in d:
                    txt_file.write(" ".join(line))
                txt_file.write("\n")
                txt_file.write("---------------------------------------------------------------------------")
                txt_file.write("\n\n")
            
        plt.show()

    
    cluster_res_nums = []
    for i in cluster_index_dict.keys():
        cluster_res_nums.append(np.where(X[cluster_index_dict[i]].any(axis=0))[0])
    cluster_res_nums = np.array(cluster_res_nums)

    #offset the column indices by the minimum residue number
    cluster_res_nums = cluster_res_nums + min_res_id

    pymol.finish_launching(['pymol', '-q'])
    cmd.refresh()
    
    # load_protein=Button2
    Label2.configure(text='''your protein is loading ,.''')
    acc_id = str(Entry1.get())
    api_request=requests.get("https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/"+acc_id)
    # pdb_dow="https://files.rcsb.org/download/"
    api=json.loads(api_request.content)
    api=api[acc_id][0]['pdb_id']
    load_protein="https://files.rcsb.org/download/"+api+".pdb"

    cmd.load(load_protein, acc_id)
    
    #Hide all representations, then show the surface representation in white colour
    cmd.hide()
    cmd.show("surface")
    cmd.color("white")

    #Make a list of desired cluster numbers to be displayed and corresponding colours
    cltsr=int(Entry2.get())
    cluster_numbers = list(range(cltsr))
    cluster_numbers = [int(i) for i in cluster_numbers] 
    colours = ["deepsalmon","deepblue","hotpink","chartreuse", "cyan", "darksalmon","dash","deepolive","deeppurple","firebrick"]

    for i in range(len(cluster_numbers)):
        res_list = '+'.join(map(str, cluster_res_nums[cluster_numbers[i]]))
        selection_name = "cluster" + str(cluster_numbers[i])
        
        #Make selection using residue list and name it
        cmd.select(selection_name, "resi " + res_list)
        
        #Colour the selection
        cmd.color(colours[i], selection_name)
        
        #Hide the pink dots by disabling selection - but selection can be continued to be manipulated
        cmd.disable(selection_name)

    #Colour the common residues with gray
    def find_common_residues():
        common = cluster_res_nums[cluster_numbers[0]]
        for i in range(1, len(cluster_numbers)):
            common = np.intersect1d(common, cluster_res_nums[cluster_numbers[i]])  
        return common
            
    common_res = find_common_residues()
    res_list = '+'.join(map(str, common_res))
    selection_name = "common"

    #Make selection using residue list and name it
    cmd.select(selection_name, "resi " + res_list)

    #Colour the selection
    cmd.color("gray60", selection_name)
    
    #Hide the pink dots by disabling selection - but selection can be continued to be manipulated
    cmd.disable(selection_name)
    Label2.configure(text=''' ''')
    Button3.configure(text='''Try again''')
    # plot_clusters()
    # New feature-----------------------------------------
    
    # Compute optimal number of clusters to be used
    # Define a list with the different k-values to be compared
    # The maximum number of clusters possible is n_samples - 1. Take the minimum of this and the fixed upper range desired.
    n_samples = X.shape[0]
    range_n_clusters = np.arange(2, min(n_samples, 16))
    # Call the function that calculates silhouette scores for different values of k and returns the optimal number of clusters
    optimal_k = compute_optimal_k(range_n_clusters)
    # # Plot the Elbow Coefficient graph used to visually determine the optimal number of clusters
    # # plot_elbow_coefficient_graph(range_n_clusters)
    # # plot_elbow_coefficient_graph()
    print("Optimal number of clusters: ", optimal_k)
    text_label_2.configure(text='''Optimal number of clusters: '''+str(optimal_k),foreground="green")
    # text_label_2.configure(text="Optimal number of clusters:  "+str(optimal_k))
    plot_clusters()
    # # --------------------------------------------------


def check_val():
    if not len(Entry1.get()) == 0:
        uniprot_id_val = re.compile(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})')
        mo = uniprot_id_val.search(Entry1.get())
        if mo is not None:
            validation()
        else:
            text_label_1.configure(text="Please check Uniport ID") 
    else:
        text_label_1.configure(text="Please enter Uniport ID")   
            
def validation(*args, **kw):
    def get_protein_info(acc_id):
        
        url = api_url + acc_id
        time.sleep(0.01)
        response = requests.get(url)
        if response.status_code == 200:
            text_label_1.configure(text='''Valid uniprot ID''',foreground="green")
            pdb_best(acc_id)
            
            return json.loads(response.content.decode('utf-8'))  
        else:
            return None
    
    acc_id = Entry1.get()
    def pdb_best(prtn_acc):
        try:
            api_request=requests.get("https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/"+prtn_acc)
            api=json.loads(api_request.content)
        except Exception as e:
            api="error..."
            message_1.configure(text=api)

        a=[]
        for x in range(5):
            a.append("PDB "+str(x)+" : "+api[acc_id][x]['pdb_id'])
        Labelframe2.configure(text='''Best PDB structure for given Accession ID''')
        message_1.configure(text=a)
   
    json_response = get_protein_info(acc_id)
    if(json_response == None):
        text_label_1.configure(text='''Invalid uniprot ID''')
        sys.exit()

# -------------------GUI--------------------------------------------------

def mytkdialog(parent):
    global Entry1,text_label_1,Entry2,message_1,Label2,Button3,Button2,Checkbutton2,Labelframe2,cluster,message_2,text_label_2

    window=tk.Tk()
    window.geometry("800x500+351+150")
    
    _bgcolor = '#d9d9d9'  # X11 color: 'gray85'
    _fgcolor = '#000000'  # X11 color: 'black'
    _compcolor = '#d9d9d9' # X11 color: 'gray85'
    _ana1color = '#d9d9d9' # X11 color: 'gray85'
    _ana2color = '#ececec' # Closest X11 color: 'gray92'
    font10 = "-family {Segoe UI} -size 12 -weight bold"
    
    font11 = "-family {Segoe UI} -size 18 -weight bold"
    
    
    window.minsize(120, 1)
    window.maxsize(1370, 749)
    window.resizable(1, 1)
    window.title("PICv")
    window.configure(background="#fff4ea")
    window.configure(highlightbackground="#d9d9d9")
    window.configure(highlightcolor="black")

    menubar = tk.Menu(window,bg=_bgcolor,fg=_fgcolor)
    window.configure(menu = menubar)

    Image_canva = tk.Canvas(window)
    Image_canva.place(relx=-0.013, rely=0.0, relheight=0.173
                , relwidth=1.015)
    Image_canva.configure(background="#ffffff")
    # Image_canva.configure(borderwidth="2")
    Image_canva.configure(highlightbackground="#d9d9d9")
    Image_canva.configure(highlightcolor="black")
    Image_canva.configure(insertbackground="black")
    Image_canva.configure(relief="ridge")
    Image_canva.configure(selectbackground="#c4c4c4")
    Image_canva.configure(selectforeground="black")
    
    Label3 = tk.Label(window)
    Label3.place(relx=0.0, rely=0.03, height=51, width=794)
    Label3.configure(background="#fff4ea")
    Label3.configure(disabledforeground="#a3a3a3")
    Label3.configure(font=font11)
    Label3.configure(foreground="#8c52ff")
    Label3.configure(text='''Protein Interaction Clustering & Visualization''')
    
    # Button4 = tk.Label(window)
    
    # Button4.place(relx=0.0, rely=0.0, height=104, width=837)
    # Button4.configure(activebackground="#ececec")
    # Button4.configure(activeforeground="#000000")
    # Button4.configure(background="#fff4ea")
    # Button4.configure(borderwidth="0")
    # Button4.configure(disabledforeground="#a3a3a3")
    # Button4.configure(foreground="#000000")
    # Button4.configure(highlightbackground="#d9d9d9")
    # Button4.configure(highlightcolor="black")
    # photo_location = "img/picv.png"
    # global _img0
    # _img0 = tk.PhotoImage(file=photo_location)
    # Button4.configure(image=_img0)
    # Button4.configure(pady="0")

    
    accession_id = tk.Label(window)
    accession_id.place(relx=0.013, rely=0.221, height=20, width=176)
    accession_id.configure(activebackground="#f9f9f9")
    accession_id.configure(activeforeground="black")
    accession_id.configure(background="#fff4ea")
    accession_id.configure(disabledforeground="#a3a3a3")
    accession_id.configure(foreground="#0096d1")
    accession_id.configure(highlightbackground="#d9d9d9")
    accession_id.configure(highlightcolor="#000000")
    accession_id.configure(text='''Uniprot Accession ID''')

    Entry1 = tk.Entry(window)
    Entry1.place(relx=0.25, rely=0.221,height=20, relwidth=0.13)
    Entry1.configure(background="white")
    Entry1.configure(disabledforeground="#a3a3a3")
    # Entry1.configure(font="TkFixedFont")
    Entry1.configure(foreground="#000000")
    Entry1.configure(highlightbackground="#d9d9d9")
    Entry1.configure(highlightcolor="black")
    Entry1.configure(insertbackground="black")
    Entry1.configure(selectbackground="#c4c4c4")
    Entry1.configure(selectforeground="black")

    cluster = tk.Label(window)
    cluster.place(relx=0.019, rely=0.335, height=20, width=176)
    cluster.configure(activebackground="#f9f9f9")
    cluster.configure(activeforeground="black")
    cluster.configure(background="#fff4ea")
    cluster.configure(disabledforeground="#a3a3a3")
    cluster.configure(foreground="#0096d1")
    cluster.configure(highlightbackground="#d9d9d9")
    cluster.configure(highlightcolor="black")
    cluster.configure(text='''Clusters''')

    Entry2 = tk.Entry(window)
    Entry2.place(relx=0.25, rely=0.331,height=20, relwidth=0.13)
    Entry2.configure(background="white")
    Entry2.configure(disabledforeground="#a3a3a3")
    # Entry2.configure(font="TkFixedFont")
    Entry2.configure(foreground="#000000")
    Entry2.configure(highlightbackground="#d9d9d9")
    Entry2.configure(highlightcolor="black")
    Entry2.configure(insertbackground="black")
    Entry2.configure(selectbackground="#c4c4c4")
    Entry2.configure(selectforeground="black")
    
    Button1 = tk.Button(window)
    Button1.place(relx=0.413, rely=0.221, height=25, width=65)
    Button1.configure(activebackground="#ff8000")
    Button1.configure(activeforeground="#000000")
    Button1.configure(background="#a9e2f5")
    Button1.configure(borderwidth="1")
    Button1.configure(disabledforeground="#a3a3a3")
    Button1.configure(foreground="#000000")
    Button1.configure(highlightbackground="#d9d9d9")
    Button1.configure(highlightcolor="black")
    Button1.configure(pady="0")
    Button1.configure(text='''Verify''',command=check_val)
    
    text_label_1 = tk.Label(window)
    text_label_1.place(relx=0.113, rely=0.275, height=21, width=320)
    text_label_1.configure(activebackground="#f9f9f9")
    text_label_1.configure(activeforeground="black")
    text_label_1.configure(background="#fff4ea")
    text_label_1.configure(disabledforeground="#ffffff")
    text_label_1.configure(foreground="#e54b49")
    text_label_1.configure(highlightbackground="#d9d9d9")
    text_label_1.configure(highlightcolor="black")
    text_label_1.configure(justify='left')
    text_label_1.configure(text='''Ex: P61769''')
    
    text_label_2 = tk.Label(window)
    text_label_2.place(relx=0.163, rely=0.394, height=20, width=277)
    text_label_2.configure(activebackground="#f9f9f9")
    text_label_2.configure(activeforeground="black")
    text_label_2.configure(background="#fff4ea")
    text_label_2.configure(cursor="fleur")
    text_label_2.configure(disabledforeground="#a3a3a3")
    text_label_2.configure(foreground="#e54b49")
    text_label_2.configure(highlightbackground="#d9d9d9")
    text_label_2.configure(highlightcolor="black")
    text_label_2.configure(justify='left')
    text_label_2.configure(text=''' ''')

    Labelframe1 = tk.LabelFrame(window)
    Labelframe1.place(relx=0.538, rely=0.44, relheight=0.51
            , relwidth=0.415)
    Labelframe1.configure(relief='groove')
    Labelframe1.configure(foreground="#2c3d63")
    Labelframe1.configure(text='''About''')
    Labelframe1.configure(background="#a8ead5")
    Labelframe1.configure(highlightbackground="#d9d9d9")
    Labelframe1.configure(highlightcolor="black")
    
    message_2 = tk.Message(Labelframe1)
    message_2.place(relx=0.03, rely=0.122, relheight=0.624
            , relwidth=0.94, bordermode='ignore')
    message_2.configure(background="#a8ead5")
    # message_2.configure(font=font11)
    message_2.configure(foreground="#000000")
    message_2.configure(highlightbackground="#d9d9d9")
    message_2.configure(highlightcolor="black")
    message_2.configure(justify='center')
    message_2.configure(text='''This plugin is developed by Center Of Excellence Computational Genomics, R V College of Engineering, Bengaluru, India

For more information visit
www.vidyaniranjan.co.in

Cite us...''')
    message_2.configure(width=312)
    
    Labelframe2 = tk.LabelFrame(window)
    Labelframe2.place(relx=0.538, rely=0.188, relheight=0.233
            , relwidth=0.42)
    Labelframe2.configure(relief='groove')
    Labelframe2.configure(foreground="#2c3d63")
    Labelframe2.configure(labelanchor="n")
    Labelframe2.configure(text='''*****************''')
    Labelframe2.configure(background="#a8ead5")
    Labelframe2.configure(highlightbackground="#d9d9d9")
    Labelframe2.configure(highlightcolor="black")

    
    
    message_1 = tk.Message(Labelframe2)
    message_1.place(relx=0.06, rely=0.232, relheight=0.696
            , relwidth=0.887, bordermode='ignore')
    message_1.configure(background="#a8ead5")
    message_1.configure(foreground="#000000")
    message_1.configure(highlightbackground="#d9d9d9")
    message_1.configure(highlightcolor="black")
    message_1.configure(justify='center')
    message_1.configure(text='''Welcome to PICv''')
    message_1.configure(width=298)
    
    Labelframe3 = tk.LabelFrame(window)
    Labelframe3.place(relx=0.038, rely=0.46, relheight=0.269
            , relwidth=0.465)
    Labelframe3.configure(relief='groove')
    Labelframe3.configure(foreground="#2c3d63")
    Labelframe3.configure(labelanchor="n")
    Labelframe3.configure(text='''Parameters''')
    Labelframe3.configure(background="#fff4ea")
    Labelframe3.configure(highlightbackground="#d9d9d9")
    Labelframe3.configure(highlightcolor="black")
    
    value_check = tk.IntVar()
    # var=IntVar()
    Checkbutton2 = tk.Checkbutton(Labelframe3)
    Checkbutton2.place(relx=0.849, rely=0.235, relheight=0.188
            , relwidth=0.083, bordermode='ignore')
    Checkbutton2.configure(activebackground="#ececec")
    Checkbutton2.configure(activeforeground="#000000")
    Checkbutton2.configure(background="#fff4ea")
    Checkbutton2.configure(disabledforeground="#a3a3a3")
    Checkbutton2.configure(foreground="#000000")
    Checkbutton2.configure(highlightbackground="#d9d9d9")
    Checkbutton2.configure(highlightcolor="black")
    Checkbutton2.configure(justify='left')
    # Checkbutton2.configure(variable=value_check,command=disable_enable_button())
    
    

    
    Message1 = tk.Message(Labelframe3)
    Message1.place(relx=0.027, rely=0.255, relheight=0.183
                , relwidth=0.730, bordermode='ignore')
    Message1.configure(background="#fff4ea")
    Message1.configure(foreground="#0096d1")
    Message1.configure(highlightbackground="#d9d9d9")
    Message1.configure(highlightcolor="black")
    Message1.configure(text='''1) Do you want to use best PDB structure ?''')
    Message1.configure(width=270)
    
    Label1 = tk.Label(Labelframe3)
    Label1.place(relx=0.027, rely=0.416, height=23, width=254
            , bordermode='ignore')
    Label1.configure(activebackground="#f9f9f9")
    Label1.configure(activeforeground="black")
    Label1.configure(background="#fff4ea")
    Label1.configure(disabledforeground="#a3a3a3")
    Label1.configure(foreground="#e54b49")
    Label1.configure(highlightbackground="#d9d9d9")
    Label1.configure(highlightcolor="black")
    Label1.configure(text='''OR''')
    
    Message2 = tk.Message(Labelframe3)
    Message2.place(relx=0.073, rely=0.557, relheight=0.168
            , relwidth=0.618, bordermode='ignore')
    Message2.configure(background="#fff4ea")
    Message2.configure(foreground="#0096d1")
    Message2.configure(highlightbackground="#d9d9d9")
    Message2.configure(highlightcolor="black")
    Message2.configure(text='''2) Do you have a structure ?''')
    Message2.configure(width=230)
    
    Button2 = tk.Button(Labelframe3)
    Button2.place(relx=0.785, rely=0.544, height=25, width=65
                , bordermode='ignore')
    Button2.configure(activebackground="#ececec")
    Button2.configure(activeforeground="#000000")
    Button2.configure(background="#a9e2f5")
    Button2.configure(borderwidth="1")
    Button2.configure(disabledforeground="#a3a3a3")
    Button2.configure(foreground="#000000")
    Button2.configure(highlightbackground="#d9d9d9")
    Button2.configure(highlightcolor="black")
    Button2.configure(pady="0")
    Button2.configure(text='''Upload''',command=open_file)
    # uploaded_pdb = tkFileDialog.askopenfilename(title='Open PDB File', filetypes=(("pdb files", "*.pdb"),))
    
    Button3 = tk.Button(window)
    Button3.place(relx=0.038, rely=0.804, height=34, width=367)
    Button3.configure(activebackground="#ececec")
    Button3.configure(activeforeground="#000000")
    Button3.configure(background="#e9687e")
    Button3.configure(borderwidth="1")
    Button3.configure(disabledforeground="#a3a3a3")
    # Button3.configure(font=font13)
    Button3.configure(foreground="#ffffff")
    Button3.configure(highlightbackground="#d9d9d9")
    Button3.configure(highlightcolor="black")
    Button3.configure(pady="0")
    Button3.configure(text='''Visualize in Pymol''',command=check)
    
    Label2 = tk.Label(window)
    Label2.place(relx=0.038, rely=0.917, height=21, width=388)
    Label2.configure(background="#fff4ea")
    Label2.configure(disabledforeground="#a3a3a3")
    Label2.configure(foreground="#ff0000")
    Label2.configure(text='''**Follow the instructions before using it''')

    
    
        
    window.mainloop()
