# Michelle Fast, Andrew Holmes
# Bioinformatics Lab 3

import itertools
import numpy as np
import sys
import argparse


def calculate_score(src_word, matching_word):
    score = 0
    match, mismatch, gap = 5, -3, -3
    for (a, b) in zip(src_word, matching_word):
        if (a==b):
            score += match
        elif a == '.':
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

    # Assign arguments
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
    p = 0.3
    indices = (scored_words > 0).sum(axis=1).argsort()[::-1][:int(len(matching_words) * p)]

    # Store top matching words along with corresponding src words with positive score
    top_matches = []
    for index in indices:
        m_word = matching_words[index]
        s_words = np.take(src_words, [idx for idx, score in enumerate(scored_words[index]) if score > 0])
        top_matches.append((m_word, s_words))

    # Compare to database!
    VERY_BEST = {
            "seq_id": "",
            "db_sequence": "",
            "src": "",
            "score": -np.inf,
            }

    db = open("seq.txt")
    for line in db:

        # Clean up database entry
        entry_id = line.split(" ")[0]
        seq = line.split(" ")[1][:-1]

        # Loop through matching words
        for match in top_matches:

            # See if matching word is in db entry
            db_sequence_index = seq.find(match[0])
            if (db_sequence_index == -1):
                continue

            # Compare each source word with db entry
            for i, src_word in enumerate(match[1]):

                # Index of src_word in source
                query_index = src.find(src_word)

                # Define padded variables
                padded_seq = len(src) * '-' + seq + len(src) * '-'
                db_sequence_index = db_sequence_index + len(src)
                padded_src = (db_sequence_index - query_index) * ' ' + src
                padded_seq_slice = padded_seq[db_sequence_index - query_index: db_sequence_index - query_index + len(src)]

                # If no gaps
                if not gaps:

                    score = calculate_score(src, padded_seq_slice)

                    # Update overall best match
                    if score >= VERY_BEST['score']:
                        VERY_BEST['score'] = score
                        VERY_BEST['db_sequence'] = padded_seq
                        VERY_BEST['src'] = padded_src
                        VERY_BEST['seq_id'] = entry_id

                else:
                    # Score alignment
                    score = calculate_score(src, padded_seq_slice)
                    #score = calculate_score(match[1][i], match[0])

                    # Extend alignment right
                    query_right, db_sequence_right = query_index + w, db_sequence_index + w
                    run = 0
                    while(query_right < len(src) and db_sequence_right < len(seq)):

                        # Compare extension
                        if src[query_right] == padded_seq[db_sequence_right]:
                            score += 5
                            db_sequence_right += 1
                            src_word += src[query_right]
                            query_right += 1

                        else:
                            for extend in range(0,11):
                                extended_src = src[:query_right] + extend * "." + src[query_right:]
                                padded_seq_slice = padded_seq[db_sequence_index - query_index: db_sequence_index - query_index + len(extended_src)]
                                score = calculate_score(extended_src, padded_seq_slice)

                                # Update overall best match
                                if score >= VERY_BEST['score']:
                                    VERY_BEST['score'] = score
                                    VERY_BEST['db_sequence'] = padded_seq
                                    VERY_BEST['src'] = padded_src
                                    VERY_BEST['seq_id'] = entry_id

                            break


                    # Extend alignment left
                    query_left, db_sequence_left = query_index - 1, db_sequence_index - 1
                    while(query_left > 0 and db_sequence_left > 0):
                        breakpoint()

                        # Compare extension
                        try:
                            if src[query_left] == padded_seq[db_sequence_left]:
                                score += 5
                                db_sequence_left = db_sequence_left - 1
                                src_word = src[query_left] + src_word
                                query_left -= 1

                            else:
                                gapped_left = src_word
                                # gaps here?
                                break
                        except:
                            breakpoint()

                        # Update overall best match
                        if score >= VERY_BEST['score']:
                            VERY_BEST['score'] = score
                            VERY_BEST['db_sequence'] = padded_seq
                            VERY_BEST['src'] = padded_src
                            VERY_BEST['seq_id'] = entry_id


                    # Update overall best match
                    if score >= VERY_BEST['score']:
                        VERY_BEST['score'] = score
                        VERY_BEST['db_sequence'] = padded_seq
                        VERY_BEST['src'] = padded_src
                        VERY_BEST['seq_id'] = entry_id
    print_results(src, VERY_BEST)


def print_results(src, result):
    breakpoint()

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
