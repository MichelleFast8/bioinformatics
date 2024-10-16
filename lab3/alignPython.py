# Michelle Fast, Andrew Holmes
# Bioinformatics Lab 3

import itertools
import numpy as np
import sys
import argparse
from tqdm import tqdm


def calculate_score(src_word, matching_word):
    '''
    matching_word is the slice of padded_seq that aligns with src_word
    '''
    score = 0
    match, mismatch, gap = 5, -3, -3
    for (a, b) in zip(src_word, matching_word):
        if (a==b):
            score += match
        elif a == '+':
            score += gap
        else:
            score += mismatch
    return score


def main():

    # Command line arg parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', action='store_true')
    parser.add_argument('src', type=str)
    args = parser.parse_args()

    gaps = args.g
    src = args.src

    # Create words
    w = 3
    src_words = [src[i] + src[i+1] + src[i+2] for i in range(len(src)-w+1)]
    matching_words = list(itertools.product("ATCG", repeat=w))
    matching_words = ["".join(word) for word in matching_words]

    # Compare each src word to all matching words, store in matrix
    scored_words = np.zeros((len(matching_words), len(src_words))) # Rows = matching words, Cols = src words
    for j, src_word in enumerate(src_words):
        for i, matching_word in enumerate(matching_words):
            score = calculate_score(src_word, matching_word)
            scored_words[i][j] = score

    # Take only the matching_words that align best with src_words
    p = 0.2
    indices = (scored_words > 0).sum(axis=1).argsort()[::-1][:int(len(matching_words) * p)]

    # Store top matching words along with corresponding src words with positive score
    top_matches = []
    for index in indices:
        m_word = matching_words[index]
        s_words = np.take(src_words, [idx for idx, score in enumerate(scored_words[index]) if score > 0])
        top_matches.append((m_word, s_words))

    VERY_BEST = {
            "seq_id": "",
            "db_sequence": "",
            "src": "",
            "score": -np.inf,
            }

    db = open("seq.txt")
    for line in tqdm(db):

        # Clean up database entry
        entry_id = line.split(" ")[0]
        seq = line.split(" ")[1][:-1]

        # Loop through matching words
        for match in top_matches:

            # See if matching word is in db entry
            db_sequence_index = seq.find(match[0])
            if (db_sequence_index == -1):
                continue

            # Define padded sequence
            padded_seq = len(src) * '-' + seq + len(src) * '-'
            db_sequence_index = db_sequence_index + len(src)

            # Compare each source word with db entry
            for i, src_word in enumerate(match[1]):

                # Index of src_word in source
                query_index = src.find(src_word)

                # Define padded src
                padded_src = (db_sequence_index - query_index) * ' ' + src
                padded_seq_slice = padded_seq[db_sequence_index - query_index: db_sequence_index - query_index + len(src)]

                '''
                Align src with sequence where the src_word and matching word line up
                '''
                score = calculate_score(src, padded_seq_slice)

                # Update overall best match
                if score >= VERY_BEST['score']:
                    VERY_BEST['score'] = score
                    VERY_BEST['db_sequence'] = padded_seq
                    VERY_BEST['src'] = padded_src
                    VERY_BEST['seq_id'] = entry_id


                ''' Introduce Gaps '''
                if (gaps):

                    # Extend alignment right
                    src_end_match = query_index + w
                    seq_end_match = db_sequence_index + w
                    while(src_end_match < len(src) and seq_end_match < len(seq)):

                        # Compare extension: take next base if matches
                        if src[src_end_match] == padded_seq[seq_end_match]:
                            score += 5
                            src_word += src[src_end_match]
                            seq_end_match += 1
                            src_end_match += 1

                        # Compare extension: introduce gap if mismatch
                        else:

                            # Determine maximum gaps that is net positive effect
                            max_gaps = min(10, (len(src) - src_end_match) * 5 // 3)

                            for extend in range(1, max(2, max_gaps)):
                                extended_src = src[:src_end_match] + extend * "+" + src[src_end_match:] # Introduce gaps to the right of matching word
                                padded_seq_slice = padded_seq[db_sequence_index - query_index: db_sequence_index - query_index + len(extended_src)] # Grab aligned portion of sequence
                                score = calculate_score(extended_src, padded_seq_slice)

                                # Update overall best match
                                if score >= VERY_BEST['score']:
                                    VERY_BEST['score'] = score
                                    VERY_BEST['db_sequence'] = padded_seq
                                    VERY_BEST['src'] = ((db_sequence_index - query_index)*' ' + extended_src).replace('+', ' ')
                                    VERY_BEST['seq_id'] = entry_id

                            break

                    # Extend alignment left
                    src_start_match = query_index - 1
                    seq_start_match = db_sequence_index - 1

                    while(src_start_match >= 0 and seq_start_match >= 0):

                        # Compare extension
                        if src[src_start_match] == padded_seq[seq_start_match]:
                            score += 5
                            src_word = src[src_start_match] + src_word
                            seq_start_match -= 1
                            src_start_match -= 1

                        else:

                            # Determine maximum gaps that is net positive effect
                            max_gaps = min(10, (src_start_match) * 5 // 3)

                            for extend in range(1, max(2, max_gaps)):
                                extended_src = src[:src_start_match+1] + extend * "+" + src[src_start_match+1:] # Introduce gaps to left of matching word
                                padded_seq_slice = padded_seq[db_sequence_index - query_index - extend: db_sequence_index - query_index - extend + len(extended_src)] # Grab aligned portion of sequence
                                score = calculate_score(extended_src, padded_seq_slice)

                                # Update overall best match
                                if score >= VERY_BEST['score']:
                                    VERY_BEST['score'] = score
                                    VERY_BEST['db_sequence'] = padded_seq
                                    VERY_BEST['src'] = ((db_sequence_index - query_index - extend)*' ' + extended_src).replace('+', ' ')
                                    VERY_BEST['seq_id'] = entry_id
                            break

    print_results(src, VERY_BEST)


def print_results(src, result):

    # Print source word
    print("SOURCE:", src)

    # Print score
    print("Highest score:", result['score'])

    # Print line ID
    print("Best match line:", result['seq_id'])

    # Create matching bars
    matching_markers = ""
    for i, seq_char in enumerate(result['db_sequence']):
        if i < len(result['src']) and seq_char == result['src'][i]:
            matching_markers += "|"
        else:
            matching_markers += " "

    print(result['src'])
    print(matching_markers)
    print(result['db_sequence'].replace("-", " "))


if __name__ == "__main__":
    main()
