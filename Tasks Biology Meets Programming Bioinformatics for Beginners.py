#TASK 1


# First, create a string variable called ori that is equal to the Vibrio cholerae ori. Don't forget to enclose your string in quotes!

# Then, print the length of ori
ori ='ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC'
#The len() function returns the number of items (length) in an object.
x = len(ori)
print(x)


#TASK 2

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

#TASK 3
# Copy your PatternCount function from the previous step below this line
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

# Now, set Text equal to the ori of Vibrio cholerae and Pattern equal to "TGATCA"
Text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
Pattern = "TGATCA"

# Finally, print the result of calling PatternCount on Text and Pattern.

# Don't forget to use the notation print() with parentheses included!
print (PatternCount(Text, Pattern))


#TASK 4

#Exercise Break (1 point): Find the most frequent 2-mer of "GATCCAGATCCCCATAC". (You should solve this exercise by hand.)
#SOLUTION:
    
def PatternCount(Text, Pattern):

    count = 0

    for i in range(len(Text)-len(Pattern)+1):

        if Text[i:i+len(Pattern)] == Pattern:

            count = count+1

    return count

text = 'GATCCAGATCCCCATAC'

kmerlist = ['GA', 'GT', 'GC', 'GG', 'AG', 'AC', 'AT', 'AA',

            'TA', 'TG', 'TC', 'TT', 'CA', 'CG', 'CC', 'CT']

for Pattern in kmerlist:

    print(Pattern, PatternCount(text, Pattern))
    
    
#TASK 5
#To find the most frequent k-mers in a string Text, we would like to create a mapping like the example below for Text = "CGATATATCCATAG" and k = 3, where we map each k-mer appearing in  Text to its number of occurrences in  Text.  
#We call this structure a frequency map; once we have constructed a frequency map for  Text and k, we will be able to find the most frequent word  by finding the k-mer whose number of occurrences achieves a maximum.

#we can integrate the known function PatternCount into a function that creates a list of all the possible K-mer in the text. 

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def Kmercouter(Text, k):

    kmerlist = []

    for i in range(len(Text)-(k+1)):

        seq = Text[i:i + k]

        if seq not in kmerlist:

            kmerlist.append(seq)

    for kmers in kmerlist:

        print(kmers, PatternCount(Text, kmers))

Kmercouter('CGATATATCCATAG', 3)

#WHAT

n = len(Text)
for i in range(n-k+1):
    Pattern = Text[i:i+k]
    freq[Pattern] = 0
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    return freq

#TASK 6 
# Insert your completed FrequencyMap() function here.
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern] = freq[Pattern] + 1
    return freq


#I will use the following codes for the explanation. You don't need to memorize it or try to understand it right now. I'll explain them step by step.

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern] = freq[Pattern] + 1
    return freq

#First, you don't have anything but a sequence and this is important. You can intrinsically see CGA, GAT, ATA, and other 3-mers, but you have not defined 3-mers in the code, yet.

#So, you first define FrequencyMap: "def FrequencyMap(Text, k):". Then create an empty dictionary: "freq={}".

#Then you introduce the range in which you want to scan 3-mers:

n = len(Text)
    for i in range(n-k+1):

#Now, I want to emphasize it again. You haven't scanned anything yet. You just said "my range is 14-3+1 (12)". After we defined the range, we want to first define the 3-mers. We can't count 3-mers without doing this first. So:

Pattern = Text[i:i+k]

Here, k=3 and i is looping from 0 to 12.#Let's explain it more. Note that this is not part of the coding, just what the machine does with the coding to visualize why we are doing this.

Text[0:0+3] ->  Text[0:3] -> CGA

Text[1:1+3]   Text[1:4] -> GAT

Text[2:2+3]   Text[2:5] -> ATA

...

Text[12:12+3]   Text[12:15] -> TAG

#Now, these are the patterns, we have defined. Because we have just defined them, have not counted them yet, I am going to assign them the value "0":

 freq[Pattern] = 0

#Right now, we have freq={CGA=0, GAT=0, ATA=0, ..., TAG=0}. So, now the machine knows what we are looking for in the text (CGATATATCCATAG).

#Now, we tell the machine to count the 3-mers we have defined in the text. Remember, we have defined the 3-mers, but that was the end. We need to tell machine to scan the text (CGATATATCCATAG) again, this time to match the 3-mers we want:

for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern] = freq[Pattern] + 1

#Let's do it manually again to visualize. Again, note that this is not a part of the coding and it is done by the machine. Also note that the text is CGATATATCCATAG to avoid confusions.

for 0 in range(14-3+1):
            if Text[0:0+3] == Pattern:
                freq[Pattern] = freq[Pattern] + 1

#Now, what is Text[0:3]? Our text is CGATATATCCATAG and the Text[0:3] is CGA. What was the Pattern? Our first pattern (3-mer) was CGA. Is Text[0:3] == Pattern? In other terms, is CGA == CGA? The answer is yes, so we add +1 to freq[CGA]. Let's do it one more time.

for 1 in range(14-3+1):
            if Text[1:1+3] == Pattern:
                freq[Pattern] = freq[Pattern] + 1


#TASK 7
# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            pattern = key
            words.append(pattern)
    return words

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for j in range(n-k+1):
            if Pattern == Text[j:j+k]:
                freq[Pattern] += 1
    return freq

#TASK 8
# Copy your updated FrequentWords function (along with all required subroutines) below this line
def FrequentWords(Text, k):
    # your code here
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key]==m:
            words.append(key)
    return words
def FrequencyMap(Text, k):
    # your code here.
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for j in range(n-k+1):
            if Text[j:j+k]==Pattern:
                freq[Pattern]=freq[Pattern]+1
    return freq

# Now set Text equal to the Vibrio cholerae oriC and k equal to 10
Text='ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC'
k=10

# Finally, print the result of calling FrequentWords on Text and k.
print(FrequentWords(Text, k))

#TASK 9
#Write a function Reverse() that takes a string Pattern and returns a string formed by reversing all the letters of Pattern.  Recall that for strings x, y, and z, the notation z=x+y concatenates x and y together into the variable z.
# Input:  A string Pattern
# Output: The reverse of Pattern
#1
def Reverse(Pattern):
  reversed_Pattern=Pattern[::-1]
  return reversed_Pattern
#2
def reverse(Pattern):
    def split(word):
        return [char for char in word]

    x=split(Pattern)
    x.reverse()
 y=''.join(x)
 return(y)

#TASK 10
#Write a function Complement() that takes a DNA string Pattern and returns a string formed by taking the complement of each letter of Pattern (don't reverse the string yet).
# Input:  A DNA string Pattern
# Output: The complementary string of Pattern (with every nucleotide replaced by its complement).
def Complement(Pattern):
    basepairs = {"A":"T", "G":"C", "T":"A", "C":"G"}
    complement = ""
    for base in Pattern:
        complement += basepairs.get(base)
    return complement

#TASK 11
#Write a function ReverseComplement() to solve the Reverse Complement Problem, which is reproduced below. (Remember that we have already given you this code at the beginning of the module.) Then add ReverseComplement (and its needed subroutines) to Replication.py.
# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(Pattern):   
    return Complement(Reverse(Pattern))

# Copy your Reverse() function here.
def Reverse(Pattern):
    # your code here
    return Pattern[::-1]

# Copy your Complement() function here.
def Complement(Pattern):
    # your code here
    dict = {'A':'T','G':'C','T':'A','C':'G'}
    return "".join(dict[i] for i in Pattern)

#TASK 12
#Write a function PatternMatching that solves the Pattern Matching Problem. Then add PatternMatching to Replication.py.
# fill in your PatternMatching() function along with any subroutines that you need.
def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    for i in range(len(Genome)):
        if (Pattern == Genome[i:i+len(Pattern)]):
            positions += [i]
    return positions

#TASK 13
#Apply your solution to the Pattern Matching Problem to find all starting positions of "CTTGATCAT" in the Vibrio cholerae genome. (Give the positions in increasing order.)
# Copy your PatternMatching function below this line.
def PatternMatching (Pattern,Genome):
    positions = [] # output variable
    Gen= len(Genome)
    Pat= len(Pattern)
    # your code here
    for occurrence in range (Gen - Pat +1):
        if Pattern == Genome [occurrence : occurrence + Pat] :
            positions.append(occurrence)
    return positions
# The following lines will automatically read in the Vibrio cholerae genome for you and store it in a variable named v_cholerae
import sys                              # needed to read the genome
input = sys.stdin.read().splitlines()   #
v_cholerae = input[1]                   # store the genome as 'v_cholerae'


# Call PatternMatching with Pattern equal to "CTTGATCAT" and Genome equal to v_cholerae,
# and store the output as a variable called positions
positions= PatternMatching ("CTTGATCAT",v_cholerae) # we put the input data here


# print the positions variable
print (positions)

#TASK 14
#How many combined occurrences of "ATGATCAAG" and "CTTGATCAT" can you find in this region?
# Copy your PatternCount function below here
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

# On the following line, create a variable called Text that is equal to the oriC region from T petrophila
Text = "AACTCTATACCTCCTTTTTGTCGAATTTGTGTGATTTATAGAGAAAATCTTATTAACTGAAACTAAAATGGTAGGTTTGGTGGTAGGTTTTGTGTACATTTTGTAGTATCTGATTTTTAATTACATACCGTATATTGTATTAAATTGACGAACAATTGCATGGAATTGAATATATGCAAAACAAACCTACCACCAAACTCTGTATTGACCATTTTAGGACAACTTCAGGGTGGTAGGTTTCTGAAGCTCTCATCAATAGACTATTTTAGTCTTTACAAACAATATTACCGTTCAGATTCAAGATTCTACAACGCTGTTTTAATGGGCGTTGCAGAAAACTTACCACCTAAAATCCAGTATCCAAGCCGATTTCAGAGAAACCTACCACTTACCTACCACTTACCTACCACCCGGGTGGTAAGTTGCAGACATTATTAAAAACCTCATCAGAAGCTTGTTCAAAAATTTCAATACTCGAAACCTACCACCTGCGTCCCCTATTATTTACTACTACTAATAATAGCAGTATAATTGATCTGA"

# On the following line, create a variable called count_1 that is equal to the number of times
# that "ATGATCAAG" occurs in Text.
count_1 = PatternCount(Text, "ATGATCAAG")

# On the following line, create a variable called count_2 that is equal to the number of times
# that "CTTGATCAT" occurs in Text.
count_2 = PatternCount(Text, "CTTGATCAT")


# Finally, print the sum of count_1 and count_2
print (count_1 + count_2)


#TASK 15
#We can therefore define the following function that takes strings Genome and symbol as input and returns the symbol array of Genome corresponding to symbol.
#def SymbolArray(Genome, symbol):
    #array = {}
    #n = len(Genome)
    #ExtendedGenome = Genome + Genome[0:n//2]
    #for i in range(n):
        #array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    #return array
    
# Input: Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
symbol = "A"
Genome = "AAAAGGGG"
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array


# Reproduce the PatternCount function here.
def PatternCount(symbol,ExtendedGenome ):
    count = 0
    for i in range(0,len(ExtendedGenome)-len(symbol)+1):
        if ExtendedGenome[i:i+len(symbol)] == symbol:
            count = count+1
    return count

#TASK 16

# Input:  Motifs, a list of kmers (strings)
# Output: Count(Motifs), the count matrix of Motifs as a dictionary of lists
def Count(Motifs):
    count = {} # initializing the count matrix
    k = len(Motifs[0]) # length of the first kmer (but all are same length)
    for symbol in "ACGT":
        count[symbol] = [] # count matrix now has keys A, C, T, and G all with values of empty list
        for j in range(k):
            count[symbol].append(0) # count matrix now has keys A, C, G, and T all with values of a list of zeroes of length equal to the length of a kmer
    t = len(Motifs) # length of Motifs, a list of kmers (strings)
    for i in range(t): # for each kmer in Motifs
        for j in range(k): # for each element of the kmer
            symbol = Motifs[i][j] # assigns the key (symbol) to a nucleotide (ACGT) in Motifs
            
            #count[symbol] corresponds to the key of the dictionary count
            #count[symbol][j] corresponds to the position in the list assigned to the key
            count[symbol][j] += 1 # adds 1 to the position in the list assigned to the key
    return count


#TASK 17 

# Insert your Count(Motifs) function here from the last Code Challenge.

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Count(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
# Input: A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    # insert your code here
    profile1 = Count(Motifs)
    for key in "ACGT":
        for key in profile1:  
            profile[key] = [x / t for x in profile1[key]]
    return profile

#TASK 18

# Insert your Count(Motifs) function here.

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Count(Motifs):
    count = {} # initializing the count matrix
    k = len(Motifs[0]) # length of the first kmer
    for symbol in "ACGT":
        count[symbol] = [] 
        for j in range(k):
            count[symbol].append(0) 
    t = len(Motifs) # length of Motifs, a list of kmers (strings)
    for i in range(t): # for each kmer in Motifs
        for j in range(k): # for each element of the kmer
            symbol = Motifs[i][j] # assigns the key to a ACGT in Motifs
            count[symbol][j] += 1 
    return count
# j-th symbol of consensus string = symbol corresponding to a maximum element in column j of Count(Motifs)
def Consensus(Motifs):
    k = len(Motifs[0]) # set k equal to the length of Motifs[0]
    count = Count(Motifs) # set count equal to the count matrix of Motifs.
    consensus = ""  # Initialise an empty consensus string
    for j in range(k): 
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m: # range through each column of the count matrix,maximum element of column j at step j is added to m
                m = count[symbol][j]
                freq_Symbol = symbol
        consensus += freq_Symbol
    return consensus

#TASK 19

# Copy your Consensus(Motifs) function here.

# Copy your Count(Motifs) function here.

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Count(Motifs):

    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)

    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
        
#TASK 20

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    pro = 1
    for i in range(len(Text)):
        pro = pro*Profile[Text[i]][i]
    return pro
       
    return consensus
def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    for motif in Motifs:
        for index in range(len(motif)):
            if motif[index] != consensus[index]:
                count += 1
    return count

#TASK 21

# Insert your Pr(text, profile) function here from Motifs.py.

# Write your ProfileMostProbableKmer() function here.
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
# Insert your Pr(text, profile) function here from Motifs.py.
def Pr(Text, Profile):
    p=1
    k = len(Profile["A"])
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p
# Write your ProfileMostProbableKmer() function here.
def ProfileMostProbableKmer(text,k,profile):
    p=-1
    result=text[0:k]
    for i in range(len(text)-k+1):
        seq=text[i:i+k]
        pr=Pr(seq,profile)
        if pr>p:
            p=pr
            result=seq
    return result


#TASK 22

# Copy your Score(Motifs), Count(Motifs), Profile(Motifs), and Consensus(Motifs) functions here.

# Then copy your ProfileMostProbableKmer(Text, k, Profile) and Pr(Text, Profile) functions here.

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def Count(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def Consensus(Motifs):
    # insert your code here
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus 

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    for i in 'ACTG': 
        for j in range(k): 
            profile[i][j] = profile[i][j]/t   
    return profile

def Score(Motifs):
    # Insert code here
    score = 0
    k = len(Motifs[0])
    count = Count(Motifs)
    max_symbol = Consensus(Motifs)
    sum1 = 0
    for j in range(k):
        m = 0
        for symbol in "ATCG":
            if count[symbol][j] > m:
                sum1 += count[symbol][j]
    for j in range(k):
        m = 0
        for symbol in "AGTC":
            if count[symbol][j] > m:
                m = count[symbol][j]
        score += m   
    return sum1-score

def Pr(Text, Profile):
    p=1
    k = len(Profile["A"])
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p

def ProfileMostProbablePattern(text,k,profile):
    p=-1
    result=text[0:k]
    for i in range(len(text)-k+1):
        seq=text[i:i+k]
        pr=Pr(seq,profile)
        if pr>p:
            p=pr
            result=seq
    return result

def GreedyMotifSearch(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

#TASK 23

# Copy your GreedyMotifSearch function (along with all required subroutines) from Motifs.py below this line

# Copy the ten strings occurring in the hyperlinked DosR dataset below.
Dna = ["", "", "", "", "", "", "", "", "", ""]

# set t equal to the number of strings in Dna and k equal to 15

# Call GreedyMotifSearch(Dna, k, t) and store the output in a variable called Motifs

# Print the Motifs variable

# Print Score(Motifs)
def Count(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
def Consensus(Motifs):
    # insert your code here
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = Count(Motifs)
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] = profile[i][j]/t  
    return profile
def Score(Motifs):
    # Insert code here
    score = 0
    k = len(Motifs[0])
    count = Count(Motifs)
    max_symbol = Consensus(Motifs)
    sum1 = 0
    for j in range(k):
        m = 0
        for symbol in "ATCG":
            if count[symbol][j] > m:
                sum1 += count[symbol][j]
    for j in range(k):
        m = 0
        for symbol in "AGTC":
            if count[symbol][j] > m:
                m = count[symbol][j]
        score += m  
    return sum1-score
def Pr(Text, Profile):
    p=1
    k = len(Profile["A"])
    for i in range(len(Text)):
        p=p*Profile[Text[i]][i]
    return p
 
def ProfileMostProbablePattern(text,k,profile):
    p=-1
    result=text[0:k]
    for i in range(len(text)-k+1):
        seq=text[i:i+k]
        pr=Pr(seq,profile)
        if pr>p:
            p=pr
            result=seq
    return result
def GreedyMotifSearch(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]
# set t equal to the number of strings in Dna and k equal to 15
t=len(Dna)
k=15
Motifs = GreedyMotifSearch(Dna, k, t)
print(Motifs)
print(Score(Motifs))