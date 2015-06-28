#! /usr/bin/env python2.7
"""This module contains two classes which can be used to calculate optimal local or global pairwise sequence alignments using an affine gap score.

align.py
BME 205 Fall 2014, Programming Assignment #8
December 3rd, 2014
Robert Calef

This module contains two classes intended for calculating an optimal global or local 
pairwise sequence alignment using a provided substitution matrix and affine gap cost
parameters. The classes provided in this module are listed below:

local_aligner  - A class defining an object that can be used for calculating optimal
                 local alignment scores, as well as outputting the alignment in A2M format.
                 Also can be used to score a local alignment in A2M format using the
                 given parameters.

global_aligner - A class defining an object that can be used for calculating optimal
                 global alignment scores, as well as outputting the alignment in A2M 
                 format. Also can be used to score a global alignment in A2M format using
                 the given parameters.

For a much more detailed description of the two classes, and their member attributes and
functions, please see each class' respective docstring.
"""
from __future__ import print_function
import string
import sys
import urllib2



class local_aligner:
    """An object for calculating and scoring optimal local alignments with an affine gap score.

    The local_aligner class defines an object used for calculating the optimal local
    alignment of two sequences with the given set of parameters. This class also contains
    functions for scoring alignments in A2M format.

    Each instance of a local_aligner object must be initialized with a substitution matrix
    in the format used by matblas and BLAST:
        #starts with any number of comment lines, each beginning with the '#' character.
        ...
        #more comments
            A  B  C ... X
        A   4 -2  1    -4
        B  -2  4 -1    -9 
        ...
        X  -4 -9  1     4

    The key points are:
        1. The first non-comment line is the row header containing a white-space delimited 
           list of the legal characters.
        2. Each line after the row header is a row of the substitution matrix itself. These
           lines begin with a legal character Y followed by white-space delimited 
           substitution scores for the pair where the i'th score is the substitution score
           for (Y,C_i) where C_i is the i'th character in the row header.

           WARNING: The characters beginning each row must appear in the same order 
                    as in the row header.

    Sequences passed to a local_aligner object should contain only legal characters as
    defined by the provided substitution matrix.

    Each instance of the local_aligner class contains nine attributes:

        self.sub_matrix   - A dict storing the provided substitution matrix, with tuples
                            of characters as keys. self.sub_matrix[('A','D')] is the 
                            substitution score for the characters 'A' and 'D'.

        self.open_penalty - An integer specifying the open gap penalty.

        self.extend       - An integer specifying the gap extension penalty.

        self.double       - An integer specifying the double gap penalty.

        self.alphabet     - A list containing the legal characters for sequences to be
                            aligned. Obtained from the substitution matrix row header.

        self.align_matrix - A dict storing the currently cached alignment matrix. If the
                            self.align() function has not been used with this instance
                            of a local_aligner object, then self.align_matrix will be
                            an empty dict. If self.align() has been called on two sequences
                            X and Y, then self.align_matrix[(i,j)] will be a tuple of three 
                            values (M_score, Ir_score, Ic_score):
                             M_score:  An integer giving the best score of an alignment
                                       of X_l...X_i and Y_k...Y_j ending with X_i aligned
                                       to Y_j, where X_i is the i'th character of X, and 
                                       l <= i and k <= j
                             Ir_score: An integer giving the best score of an alignment
                                       of X_l...X_i and Y_k...Y_j ending with X_i aligned
                                       to a gap.
                             Ic_score: An integer giving the best score of an alignment
                                       of X_l...X_i and Y_k...Y_j ending with Y_j aligned
                                       to a gap.

                            In addition, if the align() method has been used, then the 
                            score of an optimal local alignment will be stored in 
                            self.align_matrix[(len(X),len(Y)]
                             
                            If the align() method is called with either sequence being
                            an empty sequence, then the align_matrix will be an empty dict.
        
        self.trace_matrix - A dict storing the currently cached traceback matrix. If the
                            self.align() function has not been used with this instance
                            of a local_aligner object, then self.trace_matrix will be
                            an empty dict. If self.align() has been called on two sequences
                            X and Y, then self.trace_matrix[(i,j)] will be a tuple of three 
                            strings (M_source, Ir_source, Ic_source), each of which can 
                            have one of three values, "M", "Ir", "Ic":
                             M_source:  A string specifying the source of 'M_score' in 
                                        self.align_matrix[(i,j)], can also be None to 
                                        indicate that 'M_score' is the result of beginning
                                        a new local alignment starting with X_i and Y_j 
                                        aligned.
                             Ir_source: A string specifying the source of 'Ir_score' in 
                                        self.align_matrix[(i,j)]
                             M_source:  A string specifying the source of 'Ic_score' in 
                                        self.align_matrix[(i,j)]

                            In addition, if the align() method has been used, then the 
                            index of the cell in self.align_matrix where an optimal local
                            alignment ends (i.e. the position to begin traceback from) is
                            stored as a tuple (i,j) in
                            self.align_matrix[(len(X),len(Y))]
                            Will contain None instead of a tuple if the optimal alignment 
                            was the empty alignment (also if either sequence passed to the
                            most recent call of the align() function of the calling 
                            local_aligner object was empty).


       self.row_seq       - A string storing the currently cached sequence being used
                            as the x-values of the alignment matrix (the first sequence
                            being passed to the align() function, also the master sequence
                            in A2M output). If self.align() has not been called yet with
                            this instance of a local_align object, 'row_seq' will be None.
 
       self.col_seq       - A string storing the currently cached sequence being used
                            as the y-values of the alignment matrix (the second sequence
                            being passed to the align() function, also the slave sequence
                            in A2M output). If self.align() has not been called yet with
                            this instance of a local_align object, 'col_seq' will be None.

    The local_aligner class also defines a constructor and three public methods:

        local_aligner     - An initializer method that requires an input stream 
                            containing the desired substitution matrix in the format 
                            described above. Optionally takes in integers to be used for
                            the gap open, gap extension, and double gap penalties.

        score_a2m         - Takes in two sequences in A2M format, treating the first as 
                            the master sequence and the second as the slave, and returns 
                            the score of the alignment under the parameters used when 
                            initializing the calling local_aligner object.

        align             - Takes in two sequences containing only characters present in 
                            the substitution matrix, and calculates and stores alignment 
                            and traceback matrices for these two sequences. Returns the 
                            optimal alignment score under the given parameters.

        traceback_col_seq - Returns the second sequence passed to the align() function
                            in A2M format according to the calculated alignment matrix.
                            Will print an error and exit if called before align() has been
                            called with this instance of a local_aligner object.
    
    For more information on the methods provided, please see the docstring for the 
    individual function.
    """
    def __init__(self,subst,open_penalty=12,extend=1,double=3):
        """Constuct a new local_aligner object from a substitution matrix.

        Inputs:
            subst        - An input stream containing substitution matrix data in the 
                           format described in the local_aligner class docstring.

            open_penalty - An integer specifying the desired gap open penalty.
                            Defaults to 12.

            extend       - An integer specifying the desired gap extension penalty.
                           Defaults to 1.
            
            double       - An integer specifying the desired double gap penalty.
                           Defaults to 3.
        Output:
            local_align  - A local_aligner object initialized to the parameters provided
                           as inputs.

        This method is used to create a new local_aligner object using the substitution
        scores obtained from 'subst', in addition to the other provided parameters.
        """
        #First we initialize the empty fields of a local_aligner object    
        self.open_penalty = int(open_penalty)
        self.extend = int(extend)
        self.double = int(double)
        self.align_matrix = dict()
        self.trace_matrix = dict()
        self.row_seq = None
        self.col_seq = None
        sub_matrix = dict()
        row_header = []
        #Read through all the comment lines of the substitution matrix
        #file, stopping when we reach the row header.
        for line in subst:
            if line.startswith("#"): continue
            else:
                row_header = line.split()
                break
        #Read each line of the substitution matrix, skipping the first character
        #of each row (the column header), and then storing the score for the appropriate
        #pair of characters.
        for row_num,line in enumerate(subst):
            row = line.split()
            first_char = row[0]
            for col_num,score in enumerate(row):
                if col_num == 0: continue
                second_char = row_header[col_num-1]
                sub_matrix[(first_char,second_char)] = int(score)
        #Finally, just store the substitution matrix and set of legal characters.
        self.sub_matrix = sub_matrix
        self.alphabet = row_header

    def score_a2m(self,s1,s2):
        """Given master and slave sequences in A2M format, score the alignment as a local alignment.

        Inputs:
            s1    - A string containing the master sequence to score 's2' against. This
                    sequence should contain only characters present in the substitution
                    matrix used to initialize the calling local_aligner object. 's1' will
                    be converted to upper-case before scoring.

            s2    - A string containing the slave sequence in A2M format produced from an
                    alignment of 's1' and 's2'. All upper- and lower-case characters in
                    the sequence must be present in the substitution matrix used to 
                    initialize the calling local_aligner object.
        Output:
            score - An integer specifying the score of the alignment of 's1' and 's2' 
                    under the parameters used to initialize the calling local_aligner 
                    object.

        Errors:
                  - Any illegal upper-case characters encounters in either sequence 
                    will cause an error message to be printed followed by exiting
                    via sys.exit(1). Illegal characters are any characters not contained 
                    in the substitution matrix used to initialize the calling 
                    local_aligner object.
                 
                  - If the master and slave sequences contain different numbers of 
                    alignment columns, then an error message will be printed followed 
                    by exiting via sys.exit(1).

        score_a2m is used to calculate the score of a pairwise alignment, given two 
        sequences describing that alignment in A2M format. This function uses the 
        substitution scores and gap penalties used to initialize the calling local_aligner 
        object to score the alignment. 
        """
        #Initialize variables to control the logic of scoring, 'master_pos' to indicate
        #current position in the master sequence, 'in_gap', 'in_slave_gap', and 
        #'in_alignment' to keep track of how to score gaps, and whether or not to
        #add the calculated gap penalty to the overall score. Convert 's1' to upper-case
        #for compatibility with standard substitution matrix formats.
        master_pos = 0
        in_gap = False
        in_slave_gap = False
        in_alignment = False
        gap_score = 0
        score = 0
        s1 = string.upper(s1)
        #Overall, we want to iterate over the characters of 's2', adjusting the score
        #accordingly for each character. Since a deletion or insertion in 's2' relative
        #to 's1' might not actually be part of the local alignment, we build up a gap
        #cost and only subtract it from the score if the gap is between upper-case letters
        #(i.e. the gap is part of the local alignment)
        for character in s2:
            if character == "-":
                #If the character is a dash, this is an alignment column, so increment
                #position in the master sequence. If we're in the local alignment, then
                # change current gap score according to whether we're opening a gap, 
                #extending one, or switching from a gap in one sequence to a gap in the 
                #other.
                master_pos += 1
                if not in_alignment: continue
                if in_gap: gap_score += self.extend
                else:
                    if in_slave_gap:
                        gap_score += self.double
                        in_slave_gap = False
                    else: gap_score += self.open_penalty
                    in_gap=True
            elif character.islower():
                #Similar to the dash case, except this is not an alignment column, hence
                #we don't increment the position in the master sequence.
                if not in_alignment: continue
                if in_slave_gap: gap_score += self.extend
                else:
                    if in_gap:
                        gap_score += self.double
                        in_gap = False
                    else: gap_score += self.open_penalty
                    in_slave_gap = True
            elif character.isupper():
                #Upper-case letters indicate aligned characters in 's1' and 's2', thus
                #the first upper-case letter indicates the beginning of the local
                #alignment. Also want to subtract cost for a previous gap if we hit an
                #upper-case while in the local alignment.
                if not in_alignment: in_alignment = True
                else:
                    if in_gap or in_slave_gap: 
                        score -= gap_score
                        in_gap = in_slave_gap = False
                        gap_score = 0
                master_char = s1[master_pos]
                if master_char not in self.alphabet:
                    print("ERROR: Character in master sequence not in substitution "
                          "matrix: %s\n" % (master_char), file = sys.stderr)
                    sys.exit(1)
                if character not in self.alphabet:
                    print("ERROR: Character in slave sequence not in substitution "
                          "matrix: %s\n" % (character), file = sys.stderr)
                    sys.exit(1)
                score += self.sub_matrix[(master_char,character)]
                master_pos += 1
      
        if master_pos != len(s1):
            #As master_pos is incremented everytime we see an upper-case letter or a gap
            #character, so master_pos is the number of alignment columns in the slave 
            #sequence, so if master_pos is not equal to the length of the master sequence
            #then print an error and exit due to illegal A2M file
            print("ERROR: Master sequence and slave sequence do not contain the same "
                  "number of alignment columns.\nMaster sequence: %s\nSlave sequence: %s\n"
                  "Number of alignment columns in slave: %d" % (s1,s2,master_pos))
            sys.exit(1)
        return score

    def align(self,row_seq,col_seq):
        """Calculate the optimal local alignment of two sequences, returning the alignment score.

        Inputs:
            row_seq    - A string containing the master sequence to align 'col_seq' 
                         against.
                      
            col_seq    - A string containing the slave sequence to align to 'row_seq'.
       
            WARNING: Both 'row_seq' and 'col_seq' must contain only characters present in
                     the substitution matrix used to initialize the calling local_aligner
                     object. If illegal characters are encountered, an error will be
                     printed, and sys.exit(1) will be called. All characters will be
                     converted to uppercase prior to alignment.

        Output:
            best_score - An integer giving the score of an optimal local alignment of
                         'row_seq' and 'col_seq' under the parameters provided to the
                         calling local_aligner object. 

        The align() function uses the parameters used to initialize the calling 
        local_aligner object to calculate the alignment matrices for a local alignment
        of 'row_seq' and 'col_seq' using an affine gap penalty. The alignment matrices
        will be cached in the 'self.align_matrix' attribute of the calling local_aligner
        object, that is, after calling aligner.align(s1,s2), 'aligner.align_matrix' will
        contain the alignment matrices for the sequences s1 and s2. 

        While the alignment matrix is being calculated, a traceback matrix will be 
        generated simultaneously, and stored in 'self.trace_matrix' for later use.

        This function only calculates the alignment matrices and returns the score of an
        optimal alignment. To obtain the actual sequence of the alignment in A2M format,
        use the self.traceback_col_seq() function.

        For a more rigorous description of the meaning of the elements of 
        'self.align_matrix', see the local_aligner class docstring.
        """

        #Clear the alignment matrix of any data from a previous alignment, convert
        #sequences to upper-case for compatibility with most substitution matrices and
        #initialize 'best_score' and 'best_index' to keep track of the best alignment
        #score and the index of the cell containing the best score respectively.
        #'best_index' intialized to None as we want to check if either sequence is empty,
        # and if so store the best index as None and the best score as zero to represent
        #the empty alignment.
        row_len = len(row_seq)
        col_len = len(col_seq)
        self.align_matrix.clear()
        best_score = 0
        best_index = None
        neg_inf = float("-inf")
        row_seq = string.upper(row_seq)
        col_seq = string.upper(col_seq)
        if row_len == 0 or col_len == 0:
            self.trace_matrix[(row_len,col_len)] = best_index
            self.row_seq = row_seq
            self.col_seq =  col_seq
            return best_score
        neg_inf = float("-inf")
        row_seq = string.upper(row_seq)
        col_seq = string.upper(col_seq)
        #Set up boundary conditions, setting Ir and Ic values to negative infinity so
        #we can't start on a gap, and setting M values to the characters' substitution
        #score, indicating the beginning of an alignment.
        row_first_char = row_seq[0]
        col_first_char = col_seq[0]
        for num in xrange(row_len):
            self.align_matrix[(num,0)] = (self.sub_matrix[(row_seq[num],col_first_char)]
                ,neg_inf,neg_inf)
            self.trace_matrix[(num,0)] = (None,None,None)
        for num in xrange(col_len):
            self.align_matrix[(0,num)] = (self.sub_matrix[(col_seq[num],row_first_char)],
                neg_inf,neg_inf)
            self.trace_matrix[(0,num)] = (None,None,None)
        for row in xrange(1,row_len):
            for col in xrange(1,col_len):
                #In each cell, we need to calculate the Ic, M, and Ir scores, calculated
                #from the values of the cell above, diagonally up and left, and to the left
                #respectively. For each score to calculate, we find the max of the possible
                #values, and store a corresponding string (acting as a pointer for 
                #traceback) in self.trace_matrix.
                #For the Ic score, we have three possiblities for the new score:
                #opening a gap, extending a gap, or switching to a gap in the other
                #sequence.
                M_up,Ir_up,Ic_up = self.align_matrix[(row,col-1)]
                Ic_score = M_up - self.open_penalty
                Ic_source = "M"
                dbl_gap = Ir_up - self.double
                if  dbl_gap > Ic_score:
                    Ic_score = dbl_gap
                    Ic_source = "Ir"
                extend_score = Ic_up - self.extend
                if extend_score > Ic_score:
                    Ic_score = extend_score
                    Ic_source = "Ic"
                #The score for the M matrix (corresponding to aligning X_i and Y_j) can 
                #also just be the substitution score of X_i and Y_j, indicating that a
                #new local alignment is being started. To indicate that a new alignment
                #is beginning, None is used as the traceback string, as there is no source
                #for the score. All possible values of 'M_score' are something plus
                #self.sub_matrix[(X_i,Y_j)], so we take max of the somethings, and then
                #add the substitution score.
                M_diag,Ir_diag,Ic_diag = self.align_matrix[(row-1,col-1)]
                row_char = row_seq[row]
                col_char = col_seq[col]
                M_score = 0
                M_source = None
                if M_diag > M_score:
                    M_score = M_diag
                    M_source = "M"
                if  Ir_diag > M_score:
                    M_score = Ir_diag
                    M_source = "Ir"
                if Ic_diag > M_score:
                    M_score = Ic_diag
                    M_source = "Ic"
                M_score += self.sub_matrix[(row_char,col_char)]
                #Local alignnments can only end with X_i and Y_j aligned, so just check if
                #the M score is greater than the current best scre, if so update 
                #'best_score' and 'best_index' accordingly.
                if M_score > best_score:
                    best_score = M_score
                    best_index = (row,col)
                #The Ir score is calculated much the same as the Ic score, as there also
                #only three possibilities: opening a gap, extending a gap, or switching to
                #a gap in the other sequence.
                M_left,Ir_left,Ic_left = self.align_matrix[(row-1,col)]
                Ir_score = M_left - self.open_penalty
                Ir_source = "M"
                dbl_gap = Ic_left - self.double
                if  dbl_gap > Ir_score:
                    Ir_score = dbl_gap
                    Ir_source = "Ic"
                extend_score = Ir_left - self.extend
                if extend_score > Ir_score:
                    Ir_score = extend_score
                    Ir_source = "Ir"
                self.align_matrix[(row,col)] = (M_score,Ir_score,Ic_score)
                self.trace_matrix[(row,col)] = (M_source,Ir_source,Ic_source)
        #Finally, we store the index of the cell containing the best score, or None if no
        #local alignment with a score greater than zero was found. Also cache the two
        #aligned sequences.
        self.trace_matrix[(row_len,col_len)] = best_index
        self.row_seq = row_seq
        self.col_seq =  col_seq
        return best_score

    def traceback_col_seq(self):
        """Returns the cached slave sequence in A2M format, as determined by its alignment to the cached master sequence

        Inputs:
             None
        Outputs:
             col_a2m - A string giving the slave sequence (the second sequence passed to
                       the align() function) in A2M format, describing an optimal local
                       alignment of the slave sequence to the master sequence (the first
                       sequence passed to the align() function).

             WARNING: This function should only be called after the align() function has 
                      been called from the same local_aligner object, otherwise there's no
                      cached alignment to traceback on. If this function is called before
                      the align() function has been called, an error will be printed and
                      this function will call sys.exit(1) to indicate faulty use of this 
                      class.

        traceback_col_seq() is used to perform the traceback step of the alignment
        algorithm, returning the slave sequence in A2M format, providing a description of
        an optimal local aligment of the slave sequence to the master sequence.
        """
        #Check if an alignment has been calculated, if not, print an error and exit
        if self.row_seq is None or self.col_seq is None:
            print("ERROR: Cannot traceback alignment before calling the align() function.",
                file = sys.stderr)
            sys.exit(1)
        #Initialize 'col_a2m' to store the slave sequence in A2M format as we build it up
        #from the end. Also get the index of the cell containing the end of the optimal
        #local alignment, as a local alignment can only end with X_i aligned to Y_j, we
        #start with "M" as the current source. If the best index is None, then the best
        #alignment was the empty alignment, so the current indices are set to (-1,-1) to
        #indicate no characters.
        col_a2m = ""
        row_len = len(self.row_seq)
        col_len = len(self.col_seq)
        best_index = self.trace_matrix[(row_len,col_len)]
        if best_index is None: curr_col = curr_row = -1
        else: curr_row,curr_col = best_index
        curr_source = "M"
        #As the local alignment can end in the middle of either sequence, we need to
        #add the appropriate lower-case letters (corresponding to insertions in the
        #slave sequence) and the appropriate number of gaps (corresponding to deletions
        #in the slave sequence) before beginning the traceback. If current indices are
        #(-1,-1) then the whole slave sequence in lower-case is concatenated to a number
        #of gap characters equal to the master sequence's length.
        for index in reversed(xrange(curr_col+1,col_len)):
            col_a2m += string.lower(self.col_seq[index])
        col_a2m += ("-" * (row_len - curr_row - 1))
        while(True):
            #At each step of traceback, we add the appropriate character for the source
            #of the current score: a gap, lower-case character of the slave sequence, or
            #upper-case character of the slave sequence if the current score came from the
            #cell to the left, cell above, or the cell diagonally above and to the left
            #respectively. If the current row or column is less than zero, then we've
            #reached the end of one of the sequences, so break out of traceback.
            if curr_row < 0 or curr_col < 0:
                break
            M_source, Ir_source, Ic_source = self.trace_matrix[(curr_row,curr_col)]
            seq_char = self.col_seq[curr_col]
            if curr_source is None:
                #curr_source being None indicates that the local alignment began here, so
                #break out of the traceback loop.
                break
            elif curr_source == "M":
                next_char = seq_char
                next_source = M_source
                curr_row -=1
                curr_col -=1
            elif curr_source == "Ir":
                next_char = "-"
                next_source = Ir_source
                curr_row -= 1
            elif curr_source == "Ic":
                next_char = string.lower(seq_char)
                next_source = Ic_source
                curr_col -= 1
            col_a2m += next_char
            curr_source = next_source
        #After we exit the traceback, need to add in any beginning gaps or insertions just
        #as we added end gaps and insertions before the traceback. Return the reverse of
        #the A2M sequence being built up, since it was being built from the end.
        for index in reversed(xrange(0,curr_col+1)):
            col_a2m += string.lower(self.col_seq[index])
        col_a2m += ("-" * (curr_row+1))
        return col_a2m[::-1]

class global_aligner:
    """An object for calculating and scoring optimal global alignments with an affine gap score.

    The global_aligner class defines an object used for calculating the optimal global
    alignment of two sequences with the given set of parameters. This class also contains
    functions for scoring global alignments in A2M format.

    Each instance of a global_aligner object must be initialized with a substitution matrix
    in the format used by matblas and BLAST:

        #starts with any number of comment lines, each beginning with the '#' character.
        ...
        #more comments
            A  B  C ... X
        A   4 -2  1    -4
        B  -2  4 -1    -9
        ...
        X  -4 -9  1     4

    The key points are:
        1. The first non-comment line is the row header containing a white-space delimited
           list of the legal characters.
        2. Each line after the row header is a row of the substitution matrix itself. These
           lines begin with a legal character Y followed by white-space delimited
           substitution scores for the pair where the i'th score is the substitution score
           for (Y,C_i) where C_i is the i'th character in the row header.

           WARNING: The characters beginning each row must appear in the same order
                    as in the row header.

    Sequences passed to a global_aligner object should contain only legal characters as
    defined by the provided substitution matrix.

    Each instance of the global_aligner class contains nine attributes:

        self.sub_matrix   - A dict storing the provided substitution matrix, with tuples
                            of characters as keys. self.sub_matrix[('A','D')] is the
                            substitution score for the characters 'A' and 'D'.

        self.open_penalty - An integer specifying the open gap penalty.

        self.extend       - An integer specifying the gap extension penalty.

        self.double       - An integer specifying the double gap penalty.

        self.alphabet     - A list containing the legal characters for sequences to be
                            aligned. Obtained from the substitution matrix row header.

        self.align_matrix - A dict storing the currently cached alignment matrix. If the
                            self.align() function has not been used with this instance
                            of a global_aligner object, then self.align_matrix will be
                            an empty dict. If self.align() has been called on two sequences
                            X and Y, then self.align_matrix[(i,j)] will be a tuple of three
                            values (M_score, Ir_score, Ic_score):
                             M_score:  An integer giving the best score of an alignment
                                       of X_1...X_i and Y_1...Y_j ending with X_i aligned
                                       to Y_j, where X_i is the i'th character of X.
                             Ir_score: An integer giving the best score of an alignment
                                       of X_1...X_i and Y_1...Y_j ending with X_i aligned
                                       to a gap.
                             Ic_score: An integer giving the best score of an alignment
                                       of X_1...X_i and Y_1...Y_j ending with Y_j aligned
                                       to a gap.

                            In addition, if the align() method has been used, then the
                            score of an optimal global alignment will be stored in
                            self.align_matrix[(len(X),len(Y)]


        self.trace_matrix - A dict storing the currently cached traceback matrix. If the
                            self.align() function has not been used with this instance
                            of a global_aligner object, then self.trace_matrix will be
                            an empty dict. If self.align() has been called on two sequences
                            X and Y, then self.trace_matrix[(i,j)] will be a tuple of three
                            strings (M_source, Ir_source, Ic_source), each of which can
                            have one of three values, "M", "Ir", "Ic":
                             M_source:  A string specifying the source of 'M_score' in
                                        self.align_matrix[(i,j)]
                             Ir_source: A string specifying the source of 'Ir_score' in
                                        self.align_matrix[(i,j)]
                             M_source:  A string specifying the source of 'Ic_score' in
                                        self.align_matrix[(i,j)]

                            In addition, if the align() method has been used, then the
                            source of the greatest final alignment score will be stored in
                            self.align_matrix[(len(X),len(Y))]
                            as either "M", "Ir", or "Ic" to indicate whether the optimal
                            alignment ends with characters aligned, or the final character
                            from one of the sequences aligned to a gap.
                            self.align_matrix[(len(X),len(Y))]
                            Will contain None instead of a tuple if the optimal alignment
                            was the empty alignment, or if either sequence passed to the
                            most recent call of the align() function of the calling
                            global_aligner object was empty.


       self.row_seq       - A string storing the currently cached sequence being used
                            as the x-values of the alignment matrix (the first sequence
                            being passed to the align() function, also the master sequence
                            in A2M output). If self.align() has not been called yet with
                            this instance of a global_align object, 'row_seq' will be None.

       self.col_seq       - A string storing the currently cached sequence being used
                            as the y-values of the alignment matrix (the second sequence
                            being passed to the align() function, also the slave sequence
                            in A2M output). If self.align() has not been called yet with
                            this instance of a global_align object, 'col_seq' will be None.

    The global_aligner class also defines a constructor and three public methods:

        global_aligner     - An initializer method that requires an input stream
                            containing the desired substitution matrix in the format
                            described above. Optionally takes in integers to be used for
                            the gap open, gap extension, and double gap penalties.

        score_a2m         - Takes in two sequences in A2M format, treating the first as
                            the master sequence and the second as the slave, and returns
                            the score of the alignment under the parameters used when
                            initializing the calling global_aligner object.

        align             - Takes in two sequences containing only characters present in
                            the substitution matrix, and calculates and stores alignment
                            and traceback matrices for these two sequences. Returns the
                            optimal global alignment score under the given parameters.

        traceback_col_seq - Returns the second sequence passed to the align() function
                            in A2M format according to the calculated alignment matrix.
                            Will print an error and exit if called before align() has been
                            called with this instance of a global_aligner object.

    For more information on the methods provided, please see the docstring for the
    individual function.
    """

    def __init__(self,subst,open_penalty=12,extend=1,double=3):
        """Constuct a new global_aligner object from a substitution matrix.

        Inputs:
            subst         - An input stream containing substitution matrix data in the
                            format described in the local_aligner class docstring.

            open_penalty  - An integer specifying the desired gap open penalty.
                             Defaults to 12.

            extend        - An integer specifying the desired gap extension penalty.
                            Defaults to 1.

            double        - An integer specifying the desired double gap penalty.
                            Defaults to 3.
        Output:
            global_align  - A global_aligner object initialized to the parameters provided
                           as inputs.

        This method is used to create a new global_aligner object using the substitution
        scores obtained from 'subst', in addition to the other provided parameters.
        """
        #First we initialize the empty fields of a global_aligner object
        self.open_penalty = int(open_penalty)
        self.extend = int(extend)
        self.double = int(double)
        self.align_matrix = dict()
        self.trace_matrix = dict()
        self.row_seq= None
        self.col_seq = None
        sub_matrix = dict()
        row_header = []
        #Read through all the comment lines of the substitution matrix
        #file, stopping when we reach the row header.
        for line in subst:
            if line.startswith("#"): continue
            else:
                row_header = line.split()
                break
        #Read each line of the substitution matrix, skipping the first character
        #of each row (the column header), storing the score for the appropriate
        #pair of characters
        for row_num,line in enumerate(subst):
            row = line.split()
            for col_num,score in enumerate(row):
                if col_num == 0: continue
                first_char = row_header[row_num]
                second_char = row_header[col_num-1]
                sub_matrix[(first_char,second_char)] = int(score)
        #After reading substitution matrix, just store it and the set of legal characters.
        self.sub_matrix = sub_matrix
        self.alphabet = row_header

    def score_a2m(self,s1,s2):
        """Given master and slave sequences in A2M format, score the alignment as a global alignment.

        Inputs:
            s1    - A string containing the master sequence to score 's2' against. This
                    sequence should contain only characters present in the substitution
                    matrix used to initialize the calling local_aligner object. 's1' will
                    be converted to upper-case before scoring.

            s2    - A string containing the slave sequence in A2M format produced from an
                    alignment of 's1' and 's2'. All upper- and lower-case characters in
                    the sequence must be present in the substitution matrix used to
                    initialize the calling global_aligner object.
        Output:
            score - An integer specifying the score of the alignment of 's1' and 's2'
                    under the parameters used to initialize the calling global_aligner
                    object.

        Errors:
                  - Any illegal upper-case characters encounters in either sequence
                    will cause an error message to be printed followed by exiting
                    via sys.exit(1). Illegal characters are any characters not contained
                    in the substitution matrix used to initialize the calling
                    local_aligner object.

                  - If the master and slave sequences contain different numbers of
                    alignment columns, then an error message will be printed followed
                    by exiting via sys.exit(1).

        score_a2m is used to calculate the score of a pairwise alignment, given two
        sequences describing that alignment in A2M format. This function uses the
        substitution scores and gap penalties used to initialize the calling global_aligner
        object to score the alignment.

        The alphabet of characters given in the substitution matrix will be used to
        determine legal characters in 's1' and 's2'. If any illegal characters (other than
        the gap character '-') are encountered, an error will be printed, followed by
        exiting via sys.exit(1)
        """

        #Initialize 'master_pos' to keep track of current position in the master sequence,
        #'in_gap' and 'in_slave_gap' indicate whether the previous character represented 
        #a gap in the master or slave sequence respectively.
        master_pos = 0
        in_gap = False
        in_slave_gap = False
        slave_score = 0
        s1 = string.upper(s1)
        #Overall, we want to iterate over the characters of 's2', adjusting the score
        #accordingly for each character. 
        for character in s2:
            if character == "-":
                #If the character is a dash, this is an alignment column, so increment
                #current position in the master sequence and subtract a gap penalty
                #determined by the previous character.
                master_pos += 1
                if in_gap: slave_score -= self.extend
                else:
                    if in_slave_gap:
                        slave_score -= self.double
                        in_slave_gap = False
                    else: slave_score -= self.open_penalty
                    in_gap=True
            elif character.islower():
                #Similar to the dash case, except this is not an alignment column, hence
                #we don't increment the position in the master sequence.
                if in_slave_gap: slave_score -= self.extend
                else:
                    if in_gap:
                        slave_score -= self.double
                        in_gap = False
                    else: slave_score -= self.open_penalty
                    in_slave_gap = True
            elif character.isupper():
                #Upper-case letters indicate aligned characters in 's1' and 's2', thus
                #we just add the substitution score and increment the master sequence 
                #position
                in_gap = in_slave_gap = False
                master_char = s1[master_pos]
                if master_char not in self.alphabet:
                    print("ERROR: Character in master sequence not in substitution "
                          "matrix: %s\n" % (master_char), file = sys.stderr)
                    sys.exit(1)
                if character not in self.alphabet:
                    print("ERROR: Character in slave sequence not in substitution "
                          "matrix: %s\n" % (character), file = sys.stderr)
                    sys.exit(1)
                slave_score += self.sub_matrix[(master_char,character)]
                master_pos += 1
        if master_pos != len(s1):
            #As master_pos is incremented everytime we see an upper-case letter or a gap
            #character, so master_pos is the number of alignment columns in the slave
            #sequence, so if master_pos is not equal to the length of the master sequence
            #then print an error and exit due to illegal A2M file
            print("ERROR: Master sequence and slave sequence do not contain the same "
                  "number of alignment columns.\nMaster sequence: %s\nSlave sequence: %s\n"
                  "Number of alignment columns in slave: %d" % (s1,s2,master_pos))
            sys.exit(1)
        return slave_score

    def align(self,row_seq,col_seq):
        """Calculate the optimal global alignment of two sequences, returning the alignment score.

        Inputs:
            row_seq    - A string containing the master sequence to align 'col_seq'
                         against.

            col_seq    - A string containing the slave sequence to align to 'row_seq'.

            WARNING: Both 'row_seq' and 'col_seq' must contain only characters present in
                     the substitution matrix used to initialize the calling global_aligner
                     object. If illegal characters are encountered, an error will be
                     printed, and sys.exit(1) will be called. All characters will be
                     converted to uppercase prior to alignment.

        Output:
            best_score - An integer giving the score of an optimal global alignment of
                         'row_seq' and 'col_seq' under the parameters provided to the
                         calling local_aligner object.

        The align() function uses the parameters used to initialize the calling
        global_aligner object to calculate the alignment matrices for a global alignment
        of 'row_seq' and 'col_seq' using an affine gap penalty. The alignment matrices
        will be cached in the 'self.align_matrix' attribute of the calling global_aligner
        object, that is, after calling aligner.align(s1,s2), 'aligner.align_matrix' will
        contain the alignment matrices for the sequences s1 and s2.

        While the alignment matrix is being calculated, a traceback matrix will be
        generated simultaneously, and stored in 'self.trace_matrix' for later use.

        This function only calculates the alignment matrices and returns the score of an
        optimal alignment. To obtain the actual sequence of the alignment in A2M format,
        use the self.traceback_col_seq() function.

        For a more rigorous description of the meaning of the elements of
        'self.align_matrix', see the local_aligner class docstring.
        """
        #Clear the alignment matrix of any data from a previous alignment, convert
        #sequences to upper-case for compatibility with most substitution matrices.
        row_len = len(row_seq)
        col_len = len(col_seq)
        self.align_matrix.clear()    
        neg_inf = float("-inf")
        row_len = len(row_seq)
        col_len = len(col_seq)
        row_seq = string.upper(row_seq)
        col_seq = string.upper(col_seq)
        #Set up boundary conditions, in global alignment we have a beginning "-1" row and
        #column corresponding to aligning characters of 'row_seq' or 'col_seq' to gaps.
        #In the (-1,-1) cell we set the M score to 0 and the others to negative infinity as
        #the (0,0) cell can only corrrespond to X_0 aligned to Y_0, or opening a gap in
        #either sequence. In the other (i,-1) and (-1,j) cells, only the Ir and Ic scores 
        #are above negative infinity respectively, corresponding to beginning gaps in 
        #either sequence.
        self.align_matrix[(-1,-1)] = (0,neg_inf,neg_inf)
        row=col=-1
        for row in xrange(row_len):
            self.align_matrix[(row,-1)] = (neg_inf,
                -(self.open_penalty + row * self.extend),neg_inf)
        for col in xrange(col_len):
            self.align_matrix[(-1,col)] = (neg_inf,neg_inf,
                -(self.open_penalty + col * self.extend))
        for row in xrange(row_len):
            for col in xrange(col_len):
                #In each cell, we need to calculate the Ic, M, and Ir scores, calculated
                #from the values of the cell above, diagonally up and left, and to the left
                #respectively. For each score to calculate, we find the max of the possible
                #values, and store a corresponding string (acting as a pointer for
                #traceback) in self.trace_matrix.
                #For the Ic score, we have three possiblities for the new score:
                #opening a gap, extending a gap, or switching to a gap in the other 
                #sequence.

                M_up,Ir_up,Ic_up = self.align_matrix[(row,col-1)]
                Ic_score = M_up - self.open_penalty
                Ic_source = "M"
                dbl_gap = Ir_up - self.double
                if  dbl_gap > Ic_score:
                    Ic_score = dbl_gap
                    Ic_source = "Ir"
                extend_score = Ic_up - self.extend
                if extend_score > Ic_score:
                    Ic_score = extend_score
                    Ic_source = "Ic"
                #All possible values of 'M_score' are something plus 
                #self.sub_matrix[(X_i,Y_j)], so we take max of the somethings, and then
                #add the substitution score.
                M_diag,Ir_diag,Ic_diag = self.align_matrix[(row-1,col-1)]
                row_char = row_seq[row]
                col_char = col_seq[col]
                M_score = M_diag
                M_source = "M"
                if  Ir_diag > M_score:
                    M_score = Ir_diag
                    M_source = "Ir"
                if Ic_diag > M_score:
                    M_score = Ic_diag
                    M_source = "Ic"
                M_score += self.sub_matrix[(row_char,col_char)]
                #The Ir score is calculated much the same as the Ic score, as there also
                #only three possibilities: opening a gap, extending a gap, or switching to
                #a gap in the other sequence.
                M_left,Ir_left,Ic_left = self.align_matrix[(row-1,col)]
                Ir_score = M_left - self.open_penalty
                Ir_source = "M"
                dbl_gap = Ic_left - self.double
                if  dbl_gap > Ir_score:
                    Ir_score = dbl_gap
                    Ir_source = "Ic"
                extend_score = Ir_left - self.extend
                if extend_score > Ir_score:
                    Ir_score = extend_score
                    Ir_source = "Ir"
                self.align_matrix[(row,col)] = (M_score,Ir_score,Ic_score)
                self.trace_matrix[(row,col)] = (M_source,Ir_source,Ic_source)
        #Unlike local alignment, we are not restricted to alignments ending with X_i
        #aligned to Y_j for some i and j, so we have to choose the greatest of the last 
        #three alignment scores 'final_M', 'final_Ir', and 'final_Ic', corresponding to 
        #an alignment ending with X_m aligned to Y_n, an alignment ending with X_m aligned
        #to a gap, or an alignment ending with Y_n aligned to a gap respectively.
        final_M, final_Ir, final_Ic =  self.align_matrix[(row,col)]
        align_score = final_M
        final_source = "M"
        if final_Ir > align_score:
            align_score = final_Ir
            final_source = "Ir"
        if final_Ic > align_score:
            align_score = final_Ic
            final_source = "Ic"
        #Finally, just store the best alignment score, and which of the three final
        #alignment scores it is, as well as cache 'row_seq' and 'col_seq'
        self.align_matrix[(row_len,col_len)] = align_score
        self.trace_matrix[(row_len,col_len)] = final_source
        self.row_seq = row_seq
        self.col_seq = col_seq
        return align_score

    def traceback_col_seq(self):
        """Returns the cached slave sequence in A2M format, as determined by its alignment to the cached master sequence

        Inputs:
             None
        Outputs:
             col_a2m - A string giving the slave sequence (the second sequence passed to
                       the align() function) in A2M format, describing an optimal global
                       alignment of the slave sequence to the master sequence (the first
                       sequence passed to the align() function).

             WARNING: This function should only be called after the align() function has
                      been called from the same global_aligner object, otherwise there's no
                      cached alignment to traceback on. If this function is called before
                      the align() function has been called, an error will be printed and
                      this function will call sys.exit(1) to indicate faulty use of this
                      class.

        traceback_col_seq() is used to perform the traceback step of the alignment
        algorithm, returning the slave sequence in A2M format, providing a description of
        an optimal global aligment of the slave sequence to the master sequence.
        """
        #Check if an alignment has been calculated, if not, print an error and exit
        if self.row_seq is None or self.col_seq is None:
            print("ERROR: Cannot traceback alignment before calling the align() function.",
                file = sys.stderr)
            sys.exit(1)
        #Initialize 'col_a2m' to store the slave sequence in A2M format as we build it up
        #from the end. Also get the source of the final alignment score, indicating the
        #direction in which to begin the alignment.
        col_a2m = ""
        row_len = len(self.row_seq)
        col_len = len(self.col_seq)
        curr_source = self.trace_matrix[(row_len,col_len)]
        next_source = ""
        curr_row = row_len - 1
        curr_col = col_len - 1
        while(True):
            #At each step of traceback, we add the appropriate character for the source
            #of the current score: a gap, lower-case character of the slave sequence, or
            #upper-case character of the slave sequence if the current score came from the
            #cell to the left, cell above, or the cell diagonally above and to the left
            #respectively. If the current row or column is less than zero, then we've
            #reached the end of one of the sequences, so break out of traceback.
            if curr_row < 0 or curr_col < 0:
                break
            M_source, Ir_source, Ic_source = self.trace_matrix[(curr_row,curr_col)]
            seq_char = self.col_seq[curr_col]
            if curr_source == "M":
                next_char = seq_char
                next_source = M_source
                curr_row -= 1
                curr_col -= 1
            elif curr_source == "Ir":
                next_char = "-"
                next_source = Ir_source
                curr_row -= 1
            elif curr_source == "Ic":
                next_char = string.lower(seq_char)
                next_source = Ic_source
                curr_col -= 1
            col_a2m += next_char
            curr_source = next_source
        #After we exit the traceback, need to add in any beginning gaps or insertions. 
        #Return the reverse of the A2M sequence being built up, since it was being built 
        #from the end. 
        for index in reversed(xrange(0,curr_col+1)):
            col_a2m += string.lower(self.col_seq[index])
        col_a2m += ("-" * (curr_row+1))         
        return col_a2m[::-1]
