
# Copyright 2019 Lawrence Livermore National Security. 
# Produced at the Lawrence Livermore National Laboratory.
# LLNL-CODE-781679. All Rights reserved. See file LICENSE for details.
#
# This file is part of graph-embed. For more information and source code
# availability, see github.com/LLNL/graph-embed
#
# SPDX-License-Identifier: LGPL-2.1

# Adapted from: https://plot.ly/python/3d-network-graph/, 
#               https://community.plot.ly/t/custom-button-to-modify-data-source-in-plotly-python-offline/5915, 
#               https://plot.ly/python/custom-buttons/

from numpy import *
import numpy as np
import math

#======

import plotly
#import plotly.plotly as py
import plotly.graph_objs as go
import sys

#======


random.seed(0);


from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-graph", "--graph", dest="graphpath", help="input graph FILE", metavar="FILE")
parser.add_argument("-part", "--part", dest="partpath", help="input part FILE", metavar="FILE")
#parser.add_argument("-ball", "--ball", dest="ballpath", help="input ball FILE", metavar="FILE")
parser.add_argument("-coords", "--coordinates", dest="coordspath", help="input coords FILE", metavar="FILE")
parser.add_argument("-o", "--output", dest="outputpath", help="output FILE", metavar="FILE")
# parser.add_argument("-ingredients", dest="ingpath", help="ingredients FILE", metavar="FILE")

args = parser.parse_args()

graphpath = args.graphpath
partpath = args.partpath
#ballpath = args.ballpath
coordspath = args.coordspath
outputpath = args.outputpath
# ingpath = args.ingpath


# ingfile = open(ingpath)
# ingredients = [ing for ing in ingfile.readlines()]

#ballfile = open(ballpath)
#balls = [[float(i) for i in line.split(" ")] for line in ballfile.readlines()]
#balls = [(stuff[0], stuff[1], stuff[2], stuff[3], stuff[4], stuff[5] if len(stuff) >= 3+3 else 0.0) for stuff in balls]

coordsfile = open(coordspath)
coords = [[float(i) for i in line.split(" ") if i != "\n"] for line in coordsfile.readlines()]
coords = [coord if len(coord) > 2 else [coord[0], coord[1], 0.0] for coord in coords]
#print(coords)

partfile = open(partpath)
n, K = partfile.readline().split(" ")
n = int(n)
K = int(K)
partition_sizes = [int(i) for i in partfile.readline().strip().split(" ")]

partitions = []
for i in range(K):
    partition = []
    for j in range(partition_sizes[i]):
        partition.append([int(i) for i in partfile.readline().strip().split(" ")])
    partitions.append(partition)

interpolation = [np.zeros((n, partition_sizes[0]), dtype=int)]
for size in range(K-1):
    interpolation.append(np.zeros((partition_sizes[size], partition_sizes[size + 1]), dtype=int))

for i in range(K):
    for j in range(partition_sizes[i]):
        for node in partitions[i][j]:
            interpolation[i][node][j] = 1

graphfile = open(graphpath)
edges = [(int(line.split(" ")[0]), int(line.split(" ")[1])) for line in graphfile.readlines()]

#======

def initialColors(N):
    return random.rand(N,3)


def normalize(x):
    if x<0.0:
        return 0.0
    if x>1.0:
        return 1.0
    return x

level = len(partitions)
n_aggs = partition_sizes[level-1]
part = interpolation[0]
for i in range(level - 1):
    part = part.dot(interpolation[i+1])

colors = initialColors(n_aggs)


texts = []
vertices = []
for i in range(n_aggs):
    text = ""
    aggs = []
    for j in range(n):
        if part[j][i] == 1:
            text = text + "<br>" + "test" # ingredients[j]
            aggs.append(j)
    texts.append(text)
    vertices.append(aggs)

#fullColors = []
#for i in range(partition_sizes[level - 1]):
#    for j in range(n):
#        if part[j][i] == 1:
#            fullColors.append(colors[i])

#diff = 0.01
#print(colors)
#print(partitions[0])
#for notlevel in range(K):
#    level = K - notlevel - 1
#    partition = partitions[level]
#    size = partition_sizes[level-1] if level > 0 else n
#    newColors = [(0.0, 0.0, 0.0)]*size;
#    for A in range(partition_sizes[level]):
#        old_color = colors[A]
#        for i in partition[A]:
#            newColors[i] = (normalize(old_color[0] + diff * (1.0 - 2.0 * random.random())),
#                            normalize(old_color[1] + diff * (1.0 - 2.0 * random.random())),
#                            normalize(old_color[2] + diff * (1.0 - 2.0 * random.random())))
#    all_colors[level] = colors
#    colors = newColors
#fullColors = colors


#======

axis=dict(showbackground=False, showline=False, zeroline=False, showgrid=False, showticklabels=False, title='')

layout = go.Layout(
       title='',
       width=1600, height=1200,
       showlegend=False,
       scene=dict(xaxis=dict(axis), yaxis=dict(axis), zaxis=dict(axis)),
       margin=dict(t=100),
       #hovermode='closest',
       )

#======

# bugs with numpy's linspace
def linspace (start, end, num):
    delta = 1.0 * (end - start) / (num - 1)
    seq = []
    for i in range(num):
        seq.append(start + i * delta);
    #seq.append(end)
    #seq = numpy.array(seq)
    return seq

def sphere(x, y, z, r, resolution=50):
    #just a sphere
    theta = linspace(0,2*pi,resolution)
    phi = linspace(0,pi,resolution)
    coords_x = outer(cos(theta),sin(phi))
    coords_y = outer(sin(theta),sin(phi))
    coords_z = outer(ones(resolution),cos(phi))
    return [x + r * i for i in coords_x], [y + r * i for i in coords_y], [z + r * i for i in coords_z]

#=======
    
plot_datas = []

actual_colors = []
# add nodes
for color in colors:
    actual_colors.append('rgb(' + str(int(256 * color[0])) + ',' + str(int(256 * color[1])) + ',' + str(int(256 * color[2])) + ')')

ballsize = 0
linewidth = 0
level = 0
zoomlevel = 0
if (zoomlevel == 0):
    ballsize = 3 + 2 * (1 + level)
    linewidth = 1.5 + 1.5 * (1 + level)
elif (zoomlevel == 1):
    ballsize = 3 + 2 * (1 + level)
    linewidth = 3 + 1.5 * (1 + level)
elif (zoomlevel == 2):
    ballsize = 6 + 2 * (1 + level)
    linewidth = 4.5 + 1.5 * (1 + level)
start = 0
DO_VERTICES = True
if DO_VERTICES:
    for i in range(n_aggs):
        plot_datas.append(
                go.Scatter3d(
                    x=[coords[j][0] for j in vertices[i]],  
                    y=[coords[j][1] for j in vertices[i]],  
                    z=[coords[j][2] for j in vertices[i]], 
                    visible=True,
                    mode='markers', 
                    marker=dict(
                        #symbol='circle',
                        size=ballsize,
                        opacity=1.0, 
                        color=actual_colors[i],
                        line = dict(
                            color = 'rgb(0, 0, 0)',
                            width = 1
                            ),
                        ),
                    #group=[i for i in range(n)],
                    #opacity = [0.1 for i in range(n)],
                    text=texts[i],
                    hoverinfo="text"))

    
# add edges
DO_EDGES = True
if DO_EDGES:
    Xe = []
    Ye = []
    Ze = []
    print(len(coords))
    DO_SHIFT = False
    eps = 0.0
    if DO_SHIFT:
        eps = 0.001
    for (i,j) in edges:
        Xe += [coords[i][0], coords[j][0], None]
        Ye += [coords[i][1], coords[j][1], None]
        Ze += [coords[i][2] - eps, coords[j][2] - eps, None]
    linecolor = 'rgb(75,75,75)'
    #linewidth = 4.5 #2*math.sqrt(math.sqrt(1000/E))
    lineopacity = 1 #1.0 #math.sqrt(200/E)
    plot_datas.append(go.Scatter3d(x=Xe, 
                                   y=Ye, 
                                   z=Ze,
                                   visible=True,
                                   mode='lines', 
                                   line=dict(color=linecolor, width=linewidth),
                                   hoverinfo='none', 
                                   opacity=lineopacity))

# add balls
balls = []
DO_BALLS = False
if not DO_BALLS:
    balls = []
for ball in balls:
    level, A, radius, x, y, z = ball
    level = int(level)
    A = int(A)
    d_x, d_y, d_z = sphere(x, y, z, radius)
    #print(level, A, len(all_colors[level]))
    if (len(all_colors[level]) > 0 and A != 0):
        color = all_colors[level][A]
        color = 'rgb(' + str(int(256 * color[0])) + ',' + str(int(256 * color[1])) + ',' + str(int(256 * color[2])) + ')'
        plot_datas.append(go.Surface(x=d_x, 
                                     y=d_y, 
                                     z=d_z, 
                                     visible=True,
                                     showscale=False,
                                     opacity=0.35,
                                     #surfacecolor=[colors[i+1][k] for i in range(len(d_z))],
                                     colorscale=[[0,color], [1, color]]))


begin = -1
end = 1
  
Xbe=list([begin,end,None,   begin,end,None,   begin,begin,None, end,end,None,     begin,begin,None, end,end,None,     begin,begin,None, end,end,None,   begin,end,None,   begin,end,None, begin,begin,None, end,end,None])
Ybe=list([begin,begin,None, end,end,None,     begin,end,None,   begin,end,None,   begin,begin,None, begin,begin,None, end,end,None,     end,end,None,   begin,begin,None, end,end,None,   begin,end,None,   begin,end,None])
Zbe=list([begin,begin,None, begin,begin,None, begin,begin,None, begin,begin,None, begin,end,None,   begin,end,None,   begin,end,None,   begin,end,None, end,end,None,     end,end,None,   end,end,None,     end,end,None])
borderline_opacity = 0.8
borderlines=go.Scatter3d(x=Xbe, y=Ybe, z=Zbe,
                         mode='lines',
                         line=dict(color='rgb(75,75,75)', width=5.0),
                         hoverinfo='none', opacity=borderline_opacity)

begin = 0.0
end = 10.0
Xxyz=list([begin,end,None, begin,begin,None, begin,begin,None])
Yxyz=list([begin,begin,None, begin,end,None, begin,begin,None])
Zxyz=list([begin,begin,None, begin,begin,None, begin,end,None])
xyz=go.Scatter3d(x=Xxyz, y=Yxyz, z=Zxyz,
                 mode='lines',
                 line=dict(color=['rgb(255,0,0)', 'rgb(255,0,0)', 'rgb(255,0,0)', 
                                  'rgb(0,255,0)', 'rgb(0,255,0)', 'rgb(0,255,0)', 
                                  'rgb(0,0,255)', 'rgb(0,0,255)', 'rgb(0,0,255)'], width=5.0),
                 hoverinfo='none', opacity=1.0)

#=======

plot_datas = list(plot_datas)
#plot_datas.append(borderlines)
#plot_datas.append(xyz)

DO_ANISOTROPY = True
if DO_ANISOTROPY:
    pi = 3.14159
    theta = pi / 3
    phi = pi / 6
    size = 2.0
    Xa=list([-size * math.cos(theta) * math.cos(phi), size * math.cos(theta) * math.cos(phi), None])
    Ya=list([-size * math.sin(theta) * math.cos(phi), size * math.sin(theta) * math.cos(phi), None])
    Za=list([-size * math.sin(phi), size * math.sin(phi), None])
    aniso=go.Scatter3d(x=Xa, y=Ya, z=Za,
                 mode='lines',
                 line=dict(color='rgb(0,0,0)', width=5.0),
                 hoverinfo='none', opacity=1.0)
    plot_datas.append(aniso)
    Xa=list([-size * math.sin(theta), size * math.sin(theta), None])
    Ya=list([-size * math.cos(theta), size * math.cos(theta), None])
    Za=list([0, 0, None])
    aniso=go.Scatter3d(x=Xa, y=Ya, z=Za,
                 mode='lines',
                 line=dict(color='rgb(0,0,0)', width=5.0),
                 hoverinfo='none', opacity=0.0)
    plot_datas.append(aniso)

#======

print('plotting')

fig = go.Figure(data=plot_datas, layout=layout)
plotly.offline.plot(fig, filename=outputpath)
