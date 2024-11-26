# PREPRINT is available at: https://arxiv.org/abs/2411.15955
"""Copyright 2024 Melih Şahin

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License."""



from functions import *




# first of all let us get RLIM_4(42,16)
i=4
n=42
k=16
enhanced_or_not=False #Assigning it to be True would let us use the enhanced RLIM codes proposed according to the simulation results which are equivelent to RLL codes (in this demonstration we will not be using them, as we have not used them in our paper)


code_space= give_RLIM(i,n,enhanced_or_not) # this obtains RLIM_i(n)

#print(len(give_RLIM(i,n))) #this gives the number of total codes inside RLIM_i(n)

code_space = sorted(code_space, key=bir_sayar) #POWER OPTIMIZATION CONSTRAINT:  sort based on weight

# choose block length to be 16
code_space = code_space[0:2 ** k] # get RLIM_4(42,16)
code_space = sorted(code_space, key=turn_array_into_number) # sort based on the binary value to enable binary search algorithm

how_many_data_seqeunces=200 # total_data_to_be_transmitted = k*how_many_data_seqeunces bits
how_many_data_seqeunces=2

data_to_be_encoded=[]
for h in range(how_many_data_seqeunces):
    data_to_be_encoded+=generate_random_bit_array(k)

print("our data to be encoded is as follows")
print(data_to_be_encoded)

encoded_data=[]

for h1 in range(how_many_data_seqeunces):
    encoded_data+=code_space[turn_array_into_number(data_to_be_encoded[k*h1:k*(h1+1)])]

print("encoded data is as follows")
print(encoded_data)

#Simulation:
diffusion_coef=79.4 #in μm^2/s,
num_molecules_to_emit=1000
r_r=5 # radius of fully absorbing spherical receiver in μm
r_0=10 # distance between center of transmitter and receiver in μm
signal_interval=0.2 # in seconds
counting_variance_of_the_receiver=0
channel_capacity=200

detected_molecule_count_sequence = fast_simulation_binomial(diffusion_coef, round(num_molecules_to_emit), r_r, r_0, signal_interval , encoded_data, channel_capacity, counting_variance_of_the_receiver)
print("detected count seqeunce is as follows")
print(detected_molecule_count_sequence)

#Now estimate the threshold, you can also manuualy assign a value
estimated_threshold=estimate_the_threshold(i,k,code_space,r_r,r_0,diffusion_coef,signal_interval,num_molecules_to_emit,counting_variance_of_the_receiver)
print("estimated threshold is as follows")
print(estimated_threshold)

# Now do the detection:
detected = detection(detected_molecule_count_sequence,n, estimated_threshold, i,enhanced_or_not)
print("detected seqeunce is as follows")
print(detected)

#now do error correction:
error_corrected=RLIM_error_correction(detected,n,i)
print("error corrected code is as follows")
print(error_corrected)

#now do decoding
decoded=decode(error_corrected,code_space,n,k)
print("decoded bit seqeunce is as follows")
print(decoded)
print("original bit seqeunce was as follows")
print(data_to_be_encoded)

#compute bit error rate
number_of_errenous_bits=0
for a in range(len(data_to_be_encoded)):
    if data_to_be_encoded[a]!=decoded[a]:
        number_of_errenous_bits+=1
print("number of bit errors is as follows")
print(number_of_errenous_bits)
print("bit error rate (BER) is as follows")
print(number_of_errenous_bits/k*how_many_data_seqeunces)
