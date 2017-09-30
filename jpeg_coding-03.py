DEBUG = False
import matplotlib.image as mpimg
import numpy as np
import sys
import math
from math import *
import random

import matplotlib.pyplot as plt

global table_string


k_bit_depth = 255

global c
global ct
global quantum


def init_dct(quality):
	#pi = 2*atan(1)*4#3.1415927
	global quantum
	quantum=np.zeros((8,8))
	global c
	c=np.zeros((8,8))
	global ct
	ct=np.zeros((8,8))
	   
	for i in range(0,8):
		for j in range(0,8):
			quantum[i][j]=1+((1+i+j)*quality)
			
	for j in range(0,8):
		c[0][j]=1/sqrt(8)
		ct[j][0]=c[0][j]
		
	for i in range(1,8):
		for j in range(0,8):
			c[i][j]=math.sqrt(2.0/8)*cos(pi*(2*j+1)*i/(2*8))
			#print sqrt(2/8)#*cos(pi*(2*1+1)*1/(2*8))
			ct[j][i]=c[i][j]

def forward_dct(input):
	global c
	global ct
	temp = np.zeros((8,8))
	temp1= 1.0
	output = np.zeros((8,8))
	
	#matrix multiple temp, input, ct
	for i in range(0,8):
		for j in range(0,8):
			temp[i][j]=0.0
			for k in range(0,8):
				temp[i][j]+=(int((input[i][k])-k_bit_depth/2))*ct[k][j]
	
	for i in range(0,8):
		for j in range(0,8):
			temp1=0.0
			for k in range(0,8):
				temp1+=c[i][k]*temp[k][j]
			#if(temp1>127):temp1=127
			#if(temp1<-127):temp1=-127
			output[i][j]=round(temp1)
	return output
				
def inverse_dct(input):
	temp = np.zeros((8,8))
	temp1= 1.0
	global c
	global ct
	output = np.zeros((8,8))
	
	#matrix multiple temp, input, ct
	for i in range(0,8):
		for j in range(0,8):
			temp[i][j]=0.0
			for k in range(0,8):
				temp[i][j]+=input[i][k]*c[k][j]
	
	for i in range(0,8):
		for j in range(0,8):
			temp1=0.0
			for k in range(0,8):
				temp1+=ct[i][k]*temp[k][j]
			temp1+=k_bit_depth/2
			if(temp1<0):output[i][j]=0
			elif (temp1> k_bit_depth):
				output[i][j]=k_bit_depth
			else:
				output[i][j]=round(temp1)

	return output

def mask_output(input,level):
	for i in range(0,8):
		for j in range(0,8):
			if((i+j)>level):
				input[i][j]=0

def quant_output(input):
	global quantum
	for i in range(0,8):
		for j in range(0,8):
			input[i][j]=round(input[i][j]/quantum[i][j])+128#level*int(input[i][j]/level)
			
def dequant_output(input):
	global quantum
	for i in range(0,8):
		for j in range(0,8):
			input[i][j]=(input[i][j]-128)*quantum[i][j]#level*int(input[i][j]/level)
			
			
class NodeTree(object):
    def __init__(self, left=None, right=None):
        self.left = left
        self.right = right

    def children(self):
        return (self.left, self.right)

    def nodes(self):
        return (self.left, self.right)

    def __str__(self):
        return "%s_%s" % (self.left, self.right)


## Tansverse the NodeTress in every possible way to get codings
def huffmanCodeTree(node, left=True, binString=""):
    if type(node) is str:
        return {node: binString}
    (l, r) = node.children()
    d = dict()
    d.update(huffmanCodeTree(l, True, binString + "0"))
    d.update(huffmanCodeTree(r, False, binString + "1"))
    return d




k_code_length=2
def huffman_encode(string):
	output_string=""
	dictionary_string=""
	freq = {}
	
	temp_string=string
	
	while(len(temp_string)>0):
		c=temp_string[:k_code_length]
		temp_string=temp_string[k_code_length:]
		if c in freq:
			freq[c] += 1
		else:
			freq[c] = 1
	
	#for c in string:
	#	if c in freq:
	#		freq[c] += 1
	#	else:
	#		freq[c] = 1

	#Sort the frequency table based on occurrence this will also convert the
	#dict to a list of tuples
	freq = sorted(freq.items(), key=lambda x: x[1], reverse=True)

	if DEBUG:
		print " Char | Freq "
		for key, c in freq:
			print " %4r | %d" % (key, c)

	nodes = freq

	while len(nodes) > 1:
		key1, c1 = nodes[-1]
		key2, c2 = nodes[-2]
		nodes = nodes[:-2]
		node = NodeTree(key1, key2)
		nodes.append((node, c1 + c2))
		# Re-sort the list
		nodes = sorted(nodes, key=lambda x: x[1], reverse=True)

	if DEBUG:
		print "left: %s" % nodes[0][0].nodes()[0]
		print "right: %s" % nodes[0][0].nodes()[1]

	huffmanCode = huffmanCodeTree(nodes[0][0])



	print("dictionary length: "+str(len(huffmanCode)))
	output_string+=str(tobinary(len(huffmanCode),8))
	print(tobinary(len(huffmanCode),8))
	
	print " Char | Freq  | Huffman code "
	print "-----------------------------"
	for char, frequency in freq:
		print " %-4r | %5d | %12s" % (char, frequency, huffmanCode[char])
		output_string+=tobinary(int(char, 16),8)
		output_string+=tobinary(len(huffmanCode[char]),4)
		output_string+=huffmanCode[char]
	
	print(output_string)
	print(convert_bit_string(output_string))
	
	while len(string)>0:
		left=string[:k_code_length]
		string=string[k_code_length:]
		output_string+=str(huffmanCode[left])
		
	return output_string

def tohex(val):
	h = hex((int(val,2) + (1 << 16)) % (1 << 16))
	h =h[2:]
	h=h.zfill(2)
	return h
	
	
base_64_string_data=" 0123456789abcdefghijklmnopqrstuvwxyz!#%(){}[]<>+=/*:;.,~_-@$^|`"
def val_to_base_64(val):
	return(base_64_string_data[val])

def base_64_to_val(char):
	return  base_64_string_data.index(char)
	


def convert_base_64_bit_string(string):
	#global output_file
	encode_string=""
	while(len(string)>=6):
		l=string[0:6]
		string=string[6:]

		encode_string+=val_to_base_64(int(l,2))
		#output_file.write(tohex(l))
	if(len(string)>0):
		#encode_string+=tohex(string.ljust(6, '0'))
		encode_string+=val_to_base_64(int(string.ljust(6, '0'),2))
		#output_file.write(tohex(string.ljust(8, '0')))
	return encode_string

def tobinary(val,bit_len):
	if(val>=pow(2,bit_len)):raise ValueError('to binary: val exceeds len')
	return(('{0:0'+str(bit_len)+'b}').format(val))
	
	
def convert_bit_string(string):
	#global output_file
	encode_string=""
	while(len(string)>=8):
		l=string[0:8]
		string=string[8:]
		#print(tohex(l))
		encode_string+=(tohex(l))
		#output_file.write(tohex(l))
	if(len(string)>0):
		encode_string+=tohex(string.ljust(8, '0'))
		#output_file.write(tohex(string.ljust(8, '0')))
	return encode_string
	
def rle_bit_string(string,window,run_length):
	#print(string)
	encode_string=""
	last_l=""
	
	#if the bit string length is not a match for the window size, then pad with 0'sz
	while(len(string)%window!=0):
		string+="0"
	
	while(len(string)>=window):
		l=string[0:window]
		string=string[window:]	
		encode_string+=(l)
		if(l==last_l):
			#print("match")
			repeat=True
			repeat_count=0
			while(len(string)>=window and repeat):
				if(repeat_count<pow(2,run_length)-1):
					l=string[0:window]
					if(l==last_l):
						string=string[window:]	
						repeat_count+=1
					else:
						repeat=False
						l=""
				else:
					repeat=False
					l=""
			#print(repeat_count)
			encode_string+=(tobinary(repeat_count,run_length))
		
		last_l=l
	return encode_string

def load_image(image,target):
	for i in range(0,64):
		for j in range(0,64):
			c=int(k_bit_depth*image[int(j)][int(i)][0])
			target[j][i]=c


	

def	grab_block(image,x,y,w,h):
	new_array=np.zeros((w,h))
	for i in range(0,w):
		for j in range(0,h):
			new_array[i][j]=image[i+x][j+y]
	return new_array

zigzag=(
    (0, 0),
    (0, 1), (1, 0),
    (2, 0), (1, 1), (0, 2),
    (0, 3), (1, 2), (2, 1), (3, 0),
    (4, 0), (3, 1), (2, 2), (1, 3), (0, 4),
    (0, 5), (1, 4), (2, 3), (3, 2), (4, 1), (5, 0),
    (6, 0), (5, 1), (4, 2), (3, 3), (2, 4), (1, 5), (0, 6),
    (0, 7), (1, 6), (2, 5), (3, 4), (4, 3), (5, 2), (6, 1), (7, 0),
    (7, 1), (6, 2), (5, 3), (4, 4), (3, 5), (2, 6), (1, 7),
    (2, 7), (3, 6), (4, 5), (5, 4), (6, 3), (7, 2),
    (7, 3), (6, 4), (5, 5), (4, 6), (3, 7),
    (4, 7), (5, 6), (6, 5), (7, 4),
    (7, 5), (6, 6), (5, 7),
    (6, 7), (7, 6),
    (7, 7)
)	





	

def main():
	print("main")

	random.seed(1)
	global output_file
	global block_dictionary
	global c
	
	input_image = np.full( (64,64),0)
	image_last = np.full( (128,128),0)
	image_difference = np.full( (128,128),0)
	
	input_name = sys.argv[1]
    
	
	output_filename=input_name+".txt"
	output_file = open(output_filename,'w')	
    
	file_a= input_name+".png"
	source=mpimg.imread(file_a)
	load_image(source,input_image)
	init_dct(12)
	print c
	
	
	output_image=np.zeros((64,64))
    
	output_string=""
	
	for i in range(0,8):
		for j in range(0,8):
			block = grab_block(input_image,i*8,j*8,8,8)
			block = forward_dct(block)
			quant_output(block)
			
			for z in range(0,len(zigzag)):
				a=zigzag[z][0]
				b=zigzag[z][1]
				output_string+=tobinary(int(block[a][b]),8)
			
			dequant_output(block)
			decomp = inverse_dct(block)
			output_image[i*8:i*8+8,j*8:j*8+8] = decomp

	output_string = huffman_encode(convert_bit_string(output_string))

	
	output_string = rle_bit_string(output_string,6,4)

	
	output_string = convert_base_64_bit_string(output_string)
	#output_string = convert_bit_string(output_string)
	
	output_file.write("\nbase_64\n")
	output_file.write(output_string)
	print(len(output_string)/2)
	
	
		
	fig = plt.imshow(output_image, cmap='gray', interpolation="nearest")
	plt.axis('off')
	fig.axes.get_xaxis().set_visible(False)
	fig.axes.get_yaxis().set_visible(False)
	plt.savefig("out.png",bbox_inches='tight', pad_inches = 0)
		


    
main()
