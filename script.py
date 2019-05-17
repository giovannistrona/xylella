###########
#LOAD FILES AND RUN ANALYSIS
from collections import Counter
from random import randrange,random
from igraph import Graph
from numpy import exp,arange
from numpy.random import normal
from math import sqrt


net_f = open('net_ok.csv','r')
net_ = []
w_ = []
for i in net_f:
	row = i.split(',')
	net_.append(map(int,row[:2]))
	w_.append(row[2])


net_f.close()


per_f = open('node_per.csv','r')
pers = [int(round(float(i))) for i in per_f]
per_f.close()


coord_f = open('node_coord.csv','r')
coords = [map(int,i.split(',')) for i in coord_f]
coord_f.close()

size_f = open('node_size.csv','r')
sizes = [float(i) for i in size_f]
size_f.close()

x_f = open('x_index_ok.csv','r')
x_ = [float(i.split(',')[1]) for i in x_f]
x_f.close()

x_m,x_M = min(x_),max(x_)
x_ = [(i-x_m)/(x_M-x_m) for i in x_]


w_ = map(float,w_)
w_p = [exp(-dp/100.0) for dp in w_]

dict_W = dict([[tuple(sorted((net_[i]))),w_p[i]] for i in range(len(w_p))])

g = Graph(net_,directed=False)
g.vs['xf'] = x_
g.vs['coord'] = coords
g.vs['per'] = pers
g.vs['size'] = sizes
node_number = len(g.vs)



#infection progression =(1/100)*(month^2.8)
#month_ab = [0.27,0.00,0.03,1.68,2.15,3.70,5.22,6.80,10.00,1.72,0.34,0.24]
month_ab = [0.54, 0.0, 0.06, 3.36, 4.3, 7.4, 10.44, 13.6, 20.0, 3.44, 0.68, 0.48]

for rep in range(10):
	fname = str(random())[2:]+'.csv'
	out = open('./results/'+fname,'w')
	for p_inf in arange(0.1,1.1,0.1):
		for p_ldd in arange(0,0.00055,0.00005):
			ldd_count = 0.0
			inf_count = 1.0	
			g.vs['status'] = 'S'
			g.vs['t_inf'] = 0
			infs=[randrange(node_number)]
			g.vs[infs[0]]['status']='I'     #start the infection from a single node
			I_size = g.vs[infs[0]]['size']
			month = 0
			tot_time = 0
			status=Counter(g.vs['status'])
			while len(infs)>0 and tot_time<60:
				Q = set([])
				for node in infs:
					node_x,node_y = g.vs[node]['coord']
					node_t = g.vs[node]['t_inf']					
					sss = [i.index for i in g.vs[node].neighbors() if i['status']=='S']
					len_n = float(len(g.vs[node].neighbors()))
					n_phil = int(round(g.vs[node]['per']*month_ab[month]*10))
					Q|=set(sss)
					for i in sss:
						phil = int(round((n_phil-(n_phil*p_ldd))/len_n))
						ppp = 1-((1-dict_W[tuple(sorted((node,i)))]*p_inf*g.vs['xf'][i])**(phil*exp(-14.069*exp(-0.25*node_t))))
						if random()<ppp: 
							g.vs[i]['status'] = 'I'
							I_size += g.vs[i]['size']	
							inf_count+=1.0
					found_ldd = 'no'
					jump_d = abs(normal(50000,150000))
					while found_ldd == 'no':
						ldd = randrange(node_number)
						ldd_x,ldd_y = g.vs[ldd]['coord']
						ldd_d = sqrt(((node_x-ldd_x)**2)+((node_y-ldd_y)**2))
						if ldd_d<jump_d:
							found_ldd = 'yes'
					if random()<1-((1-p_inf*g.vs['xf'][ldd])**(int(round(n_phil*p_ldd))*exp(-14.069*exp(-0.25*node_t)))):
						g.vs[ldd]['status'] = 'I'
						I_size += g.vs[ldd]['size']
						ldd_count+=1.0	
				for i in infs:
					if random()<exp(-14.069*exp(-1*g.vs[i]['t_inf'])):	#R twice I 
						g.vs[i]['status'] = 'R'
				infs=[i.index for i in g.vs.select(status='I')]
				for i in infs:
					g.vs[i]['t_inf']+=1
				month+=1
				tot_time+=1
				if month == 12:
					month = 0				
				status=Counter(g.vs['status'])
				S = status['S']/float(node_number)
				I = status['I']/float(node_number)
				R = status['R']/float(node_number)
				Q = len(Q)/float(node_number)
				Saved = 1-(I+Q+R)
				out.write(','.join(map(str,[p_inf,p_ldd,tot_time,R,I,Q,Saved,ldd_count,inf_count,I_size]))+'\n')
	out.close()

 
