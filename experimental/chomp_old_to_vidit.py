#!/usr/bin/env python
import os
import sys

def process_line (line, out_file):
    # Remove line breaks and leading white spaces.
    line = line.rstrip('\r\n')
    line = line.lstrip(' ')
    
    # Ignore boring lines.
    if len(line) == 0 or line.startswith('chain complex') or line.startswith('max dimension') :
        return
    
    if line.startswith('dimension') :
        try:
            out_file.write('-1\n')
        except:
            raise
        return
    
    if not line.startswith('boundary'):
        sys.stdout.write(line)
        raise Exception('Something is wrong 1 ... ' + token)
    
    tokens = iter(line.split())
    if tokens.next() != 'boundary':
        raise Exception('Something is wrong 2 ... ' + token)
    
    #an -> n-1
    cell_id = int( tokens.next().lstrip('a') ) - 1
    
    # = -> ignore
    if tokens.next() != '=':
        raise Exception('Something is wrong = ... ' + token)
    
    # Iterate through the linear combination.
    linear_combination = ''
    linear_combination_len = 0
    work_not_done = True
    while work_not_done:
        token = None
        linear_combination_append = ''
        
        # Get next token.
        try:
            token = tokens.next()
        except StopIteration:
            work_not_done = False
            continue
        
        # Done if token == 0.
        if token == '0':
            work_not_done = False
            continue
        # Else we assert its a sign.
        if not( token == '+' or token == '-' ):
            raise Exception('Something is wrong +- ... ' + token)
        sign = token
        
        # Get next token.
        try:
            token = tokens.next()
        except StopIteration:
            work_not_done = False
            continue
        
        # Coefficient
        linear_combination_append = token + ' ' if sign == '+' else '-' + token + ' '
        
        # Get next token.
        try:
            token = tokens.next()
        except StopIteration:
            work_not_done = False
            continue
        
        # We assert its a *.
        if not token == '*':
            raise Exception('Something is wrong * ... ' + token)
        
        # Get next token.
        try:
            token = tokens.next()
        except StopIteration:
            work_not_done = False
            continue
        
        # an -> n-1
        bdry_id = str(int( token.lstrip('a') ) - 1)
        
        # Append linear comnination .
        linear_combination_append += bdry_id
        
        linear_combination_len += 1
        linear_combination += ' ' + linear_combination_append
    
    # Write linear combination to output file.
    try:
        out_file.write( '1 ' + str(linear_combination_len) + linear_combination + '\n' )
    except:
        raise
    return

def main():
    # ToDo: Select name of input and output file.
    in_file_name  = './cache/cplx_chomp_old/cplx_g_2_m_1.hcf'
    out_file_name = './test'
    out_file = None
    
    try:
        out_file = open(out_file_name, 'w')
    except:
        raise
    
    try:
        with open(in_file_name, 'r') as in_file:
            for line in in_file:
                process_line(line, out_file)
            out_file.write('-1\n')
    except:
        raise

if __name__ == "__main__":
    main()
