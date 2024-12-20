import random

def superminhash(data, m):
    """
    SuperMinHash Algorithm Implementation.
    
    Parameters:
        data (list): The input data set (e.g., a list of integers or strings).
        m (int): The size of the MinHash signature.
    
    Returns:
        list: The SuperMinHash signature of the input data.
    """
    # Initialize signature values to infinity
    h = [float('inf')] * m
    # Auxiliary arrays for Fisher-Yates shuffle
    p = [-1] * m
    q = [-1] * m
    # Histogram and max index for optimization
    b = [0] * m
    b[m - 1] = m
    a = m - 1

    for i, element in enumerate(data):
        # Initialize pseudo-random generator with seed
        #random.seed(hash(element))
        random.seed(element)
        j = 0

        while j <= a:
            r = random.uniform(0, 1)  # Random number in [0, 1)
            k = random.randint(j, m - 1)  # Random index in [j, m - 1]

            # Initialize p[j] and p[k] if not already set
            if q[j] != i:
                q[j] = i
                p[j] = j
            if q[k] != i:
                q[k] = i
                p[k] = k

            # Swap p[j] and p[k]
            p[j], p[k] = p[k], p[j]

            # Update the signature value if needed
            if r + j < h[p[j]]:
                j_prime = min(int(h[p[j]]) if h[p[j]] != float('inf') else m, m - 1)
                h[p[j]] = r + j

                # Update histogram
                if j < j_prime:
                    b[j_prime] -= 1
                    b[j] += 1
                    while b[a] == 0:
                        a -= 1

            j += 1

    return h

def calculate_jaccard_index(signature1, signature2):
    """
    Calculate the Jaccard index from two SuperMinHash signatures.
    
    Parameters:
        signature1 (list): First SuperMinHash signature.
        signature2 (list): Second SuperMinHash signature.
    
    Returns:
        float: The estimated Jaccard index.
    """
    if len(signature1) != len(signature2):
        raise ValueError("Signatures must have the same length.")

    matches = sum(1 for x, y in zip(signature1, signature2) if x == y)
    return matches / len(signature1)

def generate_kmers(sequence, k):
    """
    Generate k-mers from a given protein sequence using a sliding window.
    
    Parameters:
        sequence (str): Protein sequence (e.g., "ACDEFGHIK").
        k (int): Length of the k-mer.
    
    Returns:
        set: A set of all k-mers generated from the sequence.
    """
    if len(sequence) < k:
        return set()  # If the sequence is shorter than k, return an empty set
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}

def calculate_real_jaccard_index(set1, set2):
    """
    Calculate the Jaccard index between two sets.
    
    Parameters:
        set1 (set): First set of k-mers.
        set2 (set): Second set of k-mers.
    
    Returns:
        float: Jaccard index between the two sets.
    """
    intersection = set1 & set2  # Intersection of the two sets
    union = set1 | set2  # Union of the two sets
    return len(intersection) / len(union) if union else 0.0

# Example usage
if __name__ == "__main__":
    # Input data

    protein1 = "MRPFFLLLLLAGALVADAFAFPIASPTGAGIARPGSPLWVHTPVAGXDAGSSAKXGGVXRKKTAVVVGAGPAGL"  # Example sequence 1
    protein2 = "MRPFFLLLLLAGALVADAFAFPIASPTGAGIARPGSPLWVHTPVAGXDAGSSAKXGGVXRKKTAVVVGAGPAGLAAALVLSRVEKEGSSGFFDRVVVLEDAPKESYDPSRAYFYNINKRGQRFTDAFGIDLTRRGLEVTEFAKRVVPADPAEVFDEEKIVRQKLSEEERKRVGTMYWIPRHELVEEIVDSIDETNERNDGTHANIELVRGXRCTHVEPTDEGLVRIVTEGKGKDETEDHLVADLCVGADGISSIVRQSLEDGRFDPAEWSNAKNPSKKFGLKKFTTPSTGLRIKGLCXQPNFAIPKGGPGXDAHEKHELENRYNYSLESATKGDTDALXLIFLPQKDPDAGAGRSVNICTMPGHDLWDTNKVRTDDGGRSAKAYFEKAFPRFDWDEIVSEDEWGLFAXTEGSRFPQCQYSPSLYVSSKPLQNGGAADGAGVVLVGDALHAFPPDLGQGVNAAFCDAMVLGESFGDAVSGALSPGDERPGSFVAGXLGSYQARNGPETRALIELARCGAPYQYNQPSRRMKLGKKLWMANVLLRLFLNKATIGLSPKPAILMMMDSRSSFRKIMRKANALTTVLWSSLLFGVASLVRSRIGV"    # Example sequence 2

    # Signature size
    k=8
    m = 200
    dataset1 = generate_kmers(protein1, k)
    dataset2 = generate_kmers(protein2, k)

    # Generate SuperMinHash signatures
    signature1 = superminhash(dataset1, m)
    signature2 = superminhash(dataset2, m)

    # Print the signatures
    print("Signature 1:", signature1)
    print("Signature 2:", signature2)

    # Calculate and print the Jaccard index
    jaccard_index = calculate_jaccard_index(signature1, signature2)
    real_jaccard_index = calculate_real_jaccard_index(dataset1, dataset2)

    print("Estimated Jaccard Index:", round(jaccard_index,3))
    print("Real Jaccard Index:", round(real_jaccard_index,3))

