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



# Definitions of Functions
import math
import random

import numpy as np

def compute_channel_coefficient(number, r_r, r_0, diff, time_in_seconds):
    if number == 1:
        return (r_r / r_0) * math.erfc((r_0 - r_r) / math.sqrt(4 * diff * time_in_seconds))
    else:
        a = (r_r / r_0) * math.erfc((r_0 - r_r) / math.sqrt(4 * diff * number * time_in_seconds))
        b = (r_r / r_0) * math.erfc((r_0 - r_r) / math.sqrt(4 * diff * (number - 1) * time_in_seconds))
        return a - b



def generate_gaussian_white_noise(variance):
    mean = 0
    std_dev = np.sqrt(variance)
    noise_value = np.random.normal(mean, std_dev)
    return round(noise_value)

def fill_an_array_with_noise(array, varaince):
    returned = []
    for a in range(len(array)):
        returned.append(array[a] + generate_gaussian_white_noise(varaince))
    return returned


def fast_simulation_binomial(diffusion_coef, M, r_r, r_0, interval_period, bit_sequence, Channel_capacity, variance_given):
    all_channel_coefficients=[]
    for j in range(Channel_capacity):
        all_channel_coefficients.append(compute_channel_coefficient(j+1,r_r,r_0,diffusion_coef,interval_period))
    molecule_count_seqeunce_to_be_returned=[]
    for a in range(len(bit_sequence)):
        min_one=min(Channel_capacity,a+1)
        our_sum=0
        for b in range(min_one):
            our_sum+=bit_sequence[a-b]*np.random.binomial(M, all_channel_coefficients[b])
        molecule_count_seqeunce_to_be_returned.append(our_sum)
    return fill_an_array_with_noise(molecule_count_seqeunce_to_be_returned,variance_given)


def give_RLIM(i,n,enhanced):
    if n < 2 * i + 1:
        raise ValueError("This implementation expects n >= 2*i + 1.")
    n=n-2
    def seqeunce_adder(seqeunces, to_be_added):
        new=[]
        for a in range(len(seqeunces)):
            j = seqeunces[a].copy()
            for b in range(len(to_be_added)):
                j.insert(b,to_be_added[b])
            new.append(j)
        return new

    Cs=[]
    for a in range(i+1):
        C=[]
        for b in range(a+2):
            sequence=[]
            for c in range(a+1):
                sequence.append(0)
            if b>0:
                sequence[-(b)]=1
            C.append(sequence)
        Cs.append(C)
    special_array=[1]
    for c in range(i):
        special_array.append(0)
    for h in range(n-i-(i-1)):
        new_1 = seqeunce_adder(Cs[-1], [0])
        new_2 = seqeunce_adder(Cs[-(1+i)], special_array)
        new = new_1 + new_2
        Cs.append(new)



    new_special=[]
    for j in range(i):
        new_special.append(0)


    full_code_space=seqeunce_adder(Cs[-1],new_special)
    to_be_removed=[]
    for j in range(n+2):
        to_be_removed.append(0)

    if not enhanced:
        full_code_space.remove(to_be_removed)

    return full_code_space

def bir_sayar(arr):
    sayı=0

    for b in range(len(arr)):
        if arr[b]==1:
            sayı+=1
    return sayı

def turn_array_into_number(array):
    number=0
    for a in range(len(array)):
        number+=array[-(a+1)]*(2**a)
    return number

def generate_random_bit_array(n):
    ret=[]
    for k in range(n):
        ret.append(random.getrandbits(1))
    return ret




def compute_expressions3(M, p1, p2, p3, p4, P0, P1, variance):
    C = M*(p1 + p2 + p3)
    D = M*(p1*(1-p1) + p2*(1-p2) + p3*(1-p3)) + variance

    A = M*(p4 + p2 + p3)
    B = M*(p4*(1-p4) + p2*(1-p2) + p3*(1-p3)) + variance

    disc = B*D * ((C - A)**2 - 2*(B - D)*math.log((math.sqrt(D)*P0) / (math.sqrt(B)*P1)))
    s = math.sqrt(disc)

    tau_plus  = (D*A - B*C + s) / (D - B)

    return tau_plus


def count_zeros_with_hat(arr,order):
    sayı = 0
    for b in range(len(arr)):
        for a in range(len(arr[b])):
            if arr[b][a]==0:
                if a>=order:
                    doğruluk=True
                    for j in range(order):
                        if arr[b][a-j-1]==1:
                            doğruluk=False
                    if doğruluk:
                        sayı+=1
    return sayı



def bir_sayar_2_lik(arr, say):
    sayı = 0
    for a in range(say):
        for b in range(len(arr[a])):
            if arr[a][b]==1:
                sayı+=1
    return sayı

def estimate_the_threshold(i,k,code_space,r_r,r_0,diff,ts,molecule_count,variance):
    channel_coefficnets=[]
    for j in range(4):
        channel_coefficnets.append(compute_channel_coefficient(1+(i+1)*j, r_r, r_0, diff, ts  ))

    P1  = bir_sayar_2_lik(code_space, 2 ** k)
    PO = count_zeros_with_hat(code_space, i)

    estimated_threshold = compute_expressions3(molecule_count, channel_coefficnets[0], channel_coefficnets[1], channel_coefficnets[2], channel_coefficnets[3], PO, P1,variance)
    return estimated_threshold

def seqeunce_assigner(sequence_of_number_of_molecules, threshold):

    bit_seqeunce = []
    for a in range(len(sequence_of_number_of_molecules)):
        if (sequence_of_number_of_molecules[a] >= threshold):
            bit_seqeunce.append(1)
        else:
            bit_seqeunce.append(0)
    return bit_seqeunce


def detection(molecule_count_seqeunce, proposed_interval, threshold_given, coding_order, enhanced):
    whole_thing = []
    sequence_of_number_of_molecules = molecule_count_seqeunce.copy()
    while (True):
        seqeunce_of_molecule_counts = sequence_of_number_of_molecules[:proposed_interval]

        threhsold =threshold_given

        returned = seqeunce_assigner(seqeunce_of_molecule_counts, threhsold)
        if not enhanced:

            are_important_ones_zeros = True

            for a in range(len(returned) - coding_order):
                if returned[a + coding_order] != 0:
                    are_important_ones_zeros = False
                    break
            if (are_important_ones_zeros):
                maxi = -1
                max_index = coding_order
                for a in range(len(seqeunce_of_molecule_counts) - coding_order):
                    if seqeunce_of_molecule_counts[a + coding_order] > maxi:
                        max_index = a + coding_order
                        maxi = seqeunce_of_molecule_counts[a + coding_order]
                returned[max_index] = 1
        whole_thing += returned
        sequence_of_number_of_molecules = sequence_of_number_of_molecules[proposed_interval:]
        if (len(sequence_of_number_of_molecules) == 0):
            return whole_thing




def RLIM_error_correction_single(arr, order):
    n = len(arr)
    out = [0] * n
    skip = order 

    for t in range(n):
        if skip > 0:
            out[t] = 0
            skip -= 1
        elif arr[t] == 1:
            out[t] = 1
            skip = order
        else:
            out[t] = 0
    return out

def RLIM_error_correction(full_seqeunce,n,i):
    now=[]
    l=round(len(full_seqeunce)/n)
    for a in range(l):
        now+=RLIM_error_correction_single(full_seqeunce[(a)*n:(a+1)*n],i)
    return now

def two_array_equality_chechker(array1,array2):
    c=0
    for a in range(len(array1)):
        if (array1[a]==array2[a]):
            c+=1
    if (c==len(array1)):
        return True
    else:
        return False

def son_biri_sıfır_yap(array):
    for a in range(len(array)):
        if array[len(array)-1-a]==1:
            array[len(array)-1-a]=0
            return array
    return array

def fast_smart_array_search(ordered_code_space, original_code_space, array, önceden_gelen=0):
    our_number=turn_array_into_number(array)
    if our_number==0:
        return 0

    uzunluk=len(ordered_code_space)
    index=int(uzunluk/2)
    ortanın_büyüklük=turn_array_into_number(ordered_code_space[index])


    if two_array_equality_chechker(ordered_code_space[index],array):
        return önceden_gelen+index
    elif uzunluk==1:
        array_new=son_biri_sıfır_yap(array)
        return fast_smart_array_search(original_code_space,original_code_space,array_new)


    elif (ortanın_büyüklük<our_number):
        return fast_smart_array_search(ordered_code_space[index+1:], original_code_space, array ,önceden_gelen+index+1)
    else:
        return fast_smart_array_search(ordered_code_space[:index], original_code_space, array ,önceden_gelen)

def turn_number_into_array(number, length):
    array=[]
    for j in range(length):
        if number>=2**(length-1-j):
            array.append(1)
            number-=2**(length-1-j)
        else:
            array.append(0)
    return array
def decode(seqeunce,code_space,n,k):
    result=[]
    how_many=len(seqeunce)//n
    for l in range(how_many):
        res = fast_smart_array_search(code_space, code_space,seqeunce[n * l:n * (l + 1)])
        result+=turn_number_into_array(res, k)
    return result

def find_best_threshold(detected_seqeunce, original_bit, M, n, i ,k,enhanced_or_not,code_space):
    error_rates=[]
    best_error=None
    for a in range(M+1):
        detected = detection(detected_seqeunce, n, a, i, enhanced_or_not)
        error_corrected = RLIM_error_correction(detected, n, i)
        decoded = decode(error_corrected, code_space, n, k)
        number_of_errenous_bits = 0
        for a in range(len(original_bit)):
            if original_bit[a] != decoded[a]:
                number_of_errenous_bits += 1
        error_rates.append(number_of_errenous_bits)
        if best_error==None or number_of_errenous_bits<best_error:
            best_error=number_of_errenous_bits
    new_best_errors=[]
    for a in range(len(error_rates)):
        if error_rates[a]==best_error:
            new_best_errors.append(a)
    new_best=int((new_best_errors[-1]+new_best_errors[0])/2)
    for j in range(len(new_best_errors)):
        if abs(new_best_errors[j]-new_best)<=1:
            return new_best
    return new_best_errors[0]

