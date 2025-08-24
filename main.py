# PREPRINT is available at: https://arxiv.org/abs/2411.15955
"""Copyright 2025 Melih Şahin

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


# With the numbers below it may take in average 30 minutes for this code to terminate
how_many_data_seqeunces=10000 # total_data_to_be_transmitted for test set in bits= k*how_many_data_seqeunces bits, the size of test set
train_data_to_be_encoded_number=2000 # total_data_to_be_transmitted for training set in bits= k*how_many_data_seqeunces bits, the size of the training set



#Simulation:
diffusion_coef=79.4 #in μm^2/s,
num_molecules_to_emit_for_uncoded=800
r_r=5 # radius of fully absorbing spherical receiver in μm
r_0=10 # distance between center of transmitter and receiver in μm
signal_interval_for_uncoded=0.2 # in seconds
counting_variance_of_the_receiver=0
channel_capacity=200

# first of all let us get RLIM_3(37,16)
i=3
n=37
k=16
enhanced_or_not=False #Assigning it to be True would let us use the enhanced RLIM codes proposed according to the simulation results which are equivelent to RLL codes (in this demonstration we will not be using them, as we have not used them in our paper)


code_space= give_RLIM(i,n,enhanced_or_not) # this obtains RLIM_i(n)
code_space = sorted(code_space, key=bir_sayar)
code_space = code_space[0:2 ** k] # get RLIM_i(n,k)
code_space = sorted(code_space, key=turn_array_into_number)




number_of_1_bits=0
for j in range(2**k):
    number_of_1_bits+=bir_sayar(code_space[j])
print(number_of_1_bits)



print("number of 1 bits of our coding method is")
print(number_of_1_bits)

print("number of 1 bits of uncoded is")
print(k*(2**(k-1)))

print("normalized molecule count is as follows")
print(round(num_molecules_to_emit_for_uncoded*(k*(2**(k-1))/number_of_1_bits)))

print("normalized signal interval is as follows")
print(signal_interval_for_uncoded*(k/n))



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

#Now estimate the threshold, you can also manuualy assign a value
estimated_threshold=estimate_the_threshold(i,k,code_space,r_r,r_0,diffusion_coef,signal_interval_for_uncoded*(k/n),round(num_molecules_to_emit_for_uncoded*(k*(2**(k-1))/number_of_1_bits)),counting_variance_of_the_receiver)
print("estimated threshold is as follows")
print(estimated_threshold)




detected_molecule_count_sequence = fast_simulation_binomial(diffusion_coef, round(num_molecules_to_emit_for_uncoded*(k*(2**(k-1))/number_of_1_bits)), r_r, r_0, signal_interval_for_uncoded*(k/n) , encoded_data, channel_capacity, counting_variance_of_the_receiver)
print("detected count seqeunce is as follows")
print(detected_molecule_count_sequence)

"below we create the train set"

train_data_to_be_encoded=[]
for h in range(train_data_to_be_encoded_number):
    train_data_to_be_encoded+=generate_random_bit_array(k)
train_encoded_data=[]
for h1 in range(train_data_to_be_encoded_number):
    train_encoded_data+=code_space[turn_array_into_number(train_data_to_be_encoded[k*h1:k*(h1+1)])]

train_detected_molecule_count_sequence = fast_simulation_binomial(diffusion_coef, round(num_molecules_to_emit_for_uncoded*(k*(2**(k-1))/number_of_1_bits)), r_r, r_0, signal_interval_for_uncoded*(k/n) , train_encoded_data, channel_capacity, counting_variance_of_the_receiver)
optimal_static_threshold=find_best_threshold(train_detected_molecule_count_sequence, train_data_to_be_encoded, round(num_molecules_to_emit_for_uncoded*(k*(2**(k-1))/number_of_1_bits)), n, i, k, enhanced_or_not, code_space)

print("optimal static threshold is as follows")
print(optimal_static_threshold)

print("now we will use optimal static threshold below ")
threshold_to_be_used=optimal_static_threshold



# Now do the detection:
detected = detection(detected_molecule_count_sequence,n, threshold_to_be_used, i,enhanced_or_not)
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
number_of_errenous_bits_optimal=0
for a in range(len(data_to_be_encoded)):
    if data_to_be_encoded[a]!=decoded[a]:
        number_of_errenous_bits_optimal+=1



print("now we will use estimated static threshold below ")
threshold_to_be_used=estimated_threshold



# Now do the detection:
detected = detection(detected_molecule_count_sequence,n, threshold_to_be_used, i,enhanced_or_not)
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
print("number of bit errors in the test set  with estimated threshold is as follows")
print(number_of_errenous_bits)
print("bit error rate (BER) in the test set  with estimated threshold is as follows")
print(number_of_errenous_bits/(k*how_many_data_seqeunces))


print("number of bit errors in the test set with optimal static threshold which is computed from the training set is as follows")
print(number_of_errenous_bits_optimal)
print("bit error rate (BER) in the test set  with optimal static threshold which is computed from the training set is as follows")
print(number_of_errenous_bits_optimal/(k*how_many_data_seqeunces))



