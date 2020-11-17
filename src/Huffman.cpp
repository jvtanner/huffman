#include "Huffman.h"
#include "Testing/HuffmanTests.h"
#include "hashmap.h"
#include "vector.h"
#include "priorityqueue.h"
#include <string>
using namespace std;

/**
 * Deallocates all nodes in a Huffman tree. We've provided this helper function
 * to you since we also use it in our test driver and figured you might want
 * to make use of it.
 */
void deleteTree(HuffmanNode* tree) {
    if (tree != nullptr) {
        deleteTree(tree->zero);
        deleteTree(tree->one);
        delete tree;
    }
}

/**
 * Constructs a Huffman coding tree for the given string, using the algorithm
 * described in lecture.
 *
 * If the input string does not contain at least two different characters,
 * this function should report an error.
 *
 * When assembling larger trees out of smaller ones, make sure to set the first
 * tree dequeued from the queue to be the zero subtree of the new tree and the
 * second tree as the one subtree.
 */
HuffmanNode* huffRecurs(PriorityQueue<HuffmanNode*>& order);

HuffmanNode* huffmanTreeFor(const string& str) {
    //Loops through the characters of the string to determine frequency
    HashMap<char, int> freqList;
    for (int i = 0; i < str.length(); i++){
        char key = str[i];
        if (freqList.containsKey(key)) {
            int frequency = freqList.get(key);
            frequency++;
            freqList.put(key, frequency);
        } else {
            freqList.put(key, 1);
        }
    }
    //makes sure that there are at least two distinct characters in the string.
    if (freqList.size() <= 1){
        error("mistake");
    }
    //Create a Huffman node for each unique character and stores them in a priority queue according to their frequency.
    PriorityQueue<HuffmanNode*> order;
    for (char key: freqList.keys()) {
        HuffmanNode* add = new HuffmanNode;
        add->one = nullptr;
        add->zero = nullptr;
        add->ch = key;
        order.enqueue(add, freqList.get(key));
    }
    return huffRecurs(order);
}
//Recursively create a Huffman tree using the priority queue that we just created.
HuffmanNode* huffRecurs(PriorityQueue<HuffmanNode*>& order){
    //base case, create a leaf node
    if (order.size() == 1) {
        HuffmanNode* result = order.dequeue();
        return result;
    //recursive case, creating the branch nodes
    } else {
        int one = order.peekPriority();
        HuffmanNode* first = order.dequeue();
        int two = order.peekPriority();
        HuffmanNode* second = order.dequeue();
        HuffmanNode* third = new HuffmanNode;
        third->ch = ' ';
        third->zero = first;
        third->one = second;
        order.enqueue(third, one+two);
        return huffRecurs(order);
    }
}


/**
 * Given a Queue<Bit> containing a compressed message and a tree that was used
 * to encode those bits, decodes the bits back to the original message.
 *
 * You can assume the input tree is not null and is not a single character;
 * those are edge cases you don't need to handle.
 *
 * You can assume the input Queue<Bit> is well-formed in that every character
 * was encoded correctly, there are no stray bits in the Queue, etc.
 */
string decodeText(Queue<Bit>& bits, HuffmanNode* tree) {
    string result;
    HuffmanNode* start = tree;
    //Follow the bits either to the left or right pointers
    while (!bits.isEmpty()) {
        Bit bit = bits.dequeue();
        if (bit == 0) {
            tree = tree->zero;
        } else {
            tree = tree->one;
        }
        //When you get the leaf, extract the value at that node
        if (tree->zero == nullptr && tree->one == nullptr) {
            result += tree->ch;
            tree = start;
        }
    }
    return result;
}

/**
 * Given a string and a Huffman encoding tree, encodes that text using the tree
 * and outputs a Queue<Bit> corresponding to the encoded representation.
 *
 * The input tree will not be null and will not consist of a single node; these
 * are edge cases you don't have to handle. The input tree will contain all
 * characters that make up the input string.
 */
void encodeRecurs(const string& str, HuffmanNode* tree, HashMap<char, Queue<Bit>>& record, Queue<Bit> trail);

//a wrapper function for the recursive strategy below
Queue<Bit> encodeText(const string& str, HuffmanNode* tree) {
    Queue<Bit> trail;
    HashMap<char, Queue<Bit>> record;
    //additional parameters are the hashmap and queue of bits
    encodeRecurs(str, tree, record, trail);
    Queue<Bit> result;
    for (char letter: str){
        Queue<Bit> name = record.get(letter);
        while (!name.isEmpty()){
            Bit bit = (name).dequeue();
            result.enqueue(bit);
        }
    }
    return result;
}

//record the bit trail involved to arrive at each leaf node, store the chars as keys and the bit trails as corresponding values
void encodeRecurs(const string& str, HuffmanNode* tree, HashMap<char, Queue<Bit>>& record, Queue<Bit> trail){
    //base case: found leaf node, add character to hashmap
    if (tree->zero == nullptr && tree->one == nullptr){
        char letter = tree->ch;
        record.put(letter, trail);
    }
    //move down the tree, adding 1 or 0 to the bit queue
    else {
        Queue<Bit> trailOne = trail;
        trailOne.enqueue(1);
        encodeRecurs(str, tree->one, record, trailOne);
        Queue<Bit> trailZero = trail;
        trailZero.enqueue(0);
        encodeRecurs(str, tree->zero, record, trailZero);
    }
}

/**
 * Encodes the given Huffman tree as a Queue<Bit> and Queue<char> in the manner
 * specified in the assignment handout.
 *
 * You can assume the input Queues are empty on entry to this function.
 *
 * You can assume that the Huffman tree provided to you is properly structured,
 * in that each internal node will have two children and only the characters in
 * the leaves matter, etc.
 */
void encodeTree(HuffmanNode* tree, Queue<Bit>& bits, Queue<char>& leaves) {
    //a leaf should enqueue a 0 for the bits, and the char for the leaves
    if (tree->zero == nullptr && tree->one == nullptr){
        bits.enqueue(0);
        char letter = tree->ch;
        leaves.enqueue(letter);
    }
    //branch should enqueue a 1 for the bits
    else {
        bits.enqueue(1);
        HuffmanNode* treeLeft = tree->zero;
        encodeTree(treeLeft, bits, leaves);
        HuffmanNode* treeRight = tree->one;
        encodeTree(treeRight, bits, leaves);
    }
}


/**
 * Decodes the given Queue<Bit> and Queue<char> into a Huffman coding tree.
 *
 * You can assume that the input Queues are structured properly in that they
 * represent a legal encoding of a tree, that there aren't stray characters
 * or bits in them, etc.
 */
HuffmanNode* decodeTree(Queue<Bit>& bits, Queue<char>& leaves) {
    //empty case returns error
    if (bits.isEmpty() || leaves.isEmpty()){
        error("nothing here");
    }
    Bit bit = bits.dequeue();
    //if the next bit is a 0, establish a leaf
    if (bit == 0){
        HuffmanNode* result = new HuffmanNode;
        result->ch = leaves.dequeue();
        result->zero = nullptr;
        result->one = nullptr;
        return result;
    }
    //if the next bit is a 1, establish a branch
    else {
        HuffmanNode* result = new HuffmanNode;
        result->zero = decodeTree(bits, leaves);
        result->one = decodeTree(bits, leaves);
        return result;
    }
}

/**
 * Compresses the given text string using Huffman coding, producing as output
 * a HuffmanResult containing the encoded tree and message.
 *
 * Your implementation of this function should report an error if there are
 * fewer than two distinct characters in the input string.
 */
HuffmanResult compress(const string& text) {
    /* TODO: Delete this comment and the next few lines, then implement this. */
    (void) text;
    return {};
}

/**
 * Decompresses the given HuffmanResult and returns the string it represents.
 *
 * Your implementation may change the file parameter however it sees fit. There
 * are no requirements about what it should look like after this function
 * returns.
 *
 * You can assume the input file is well-formed and was created by a correct
 * implementation of compress.
 */
string decompress(HuffmanResult& file) {
    /* TODO: Delete this comment and the next few lines, then implement this. */
    (void) file;
    return "";
}







/* * * * * * Test Cases Below This Point * * * * * */

ADD_TEST("Provided Test: huffmanTreeFor reports error on single type non-letter string.") {
    EXPECT_ERROR(huffmanTreeFor("!!"));    // Returns error
}


ADD_TEST("Provided Test: decodeText works on empty bit input.") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N
     *                / \
     *               M   S
     */
    HuffmanNode* tree = new HuffmanNode {
        '*',
            new HuffmanNode { 'O', nullptr, nullptr },
            new HuffmanNode { '*',
                new HuffmanNode{ '*',
                    new HuffmanNode { 'M', nullptr, nullptr },
                    new HuffmanNode { 'S', nullptr, nullptr }
                },
                new HuffmanNode{ 'N', nullptr, nullptr }
            }
    };

    /* What you get if you encode MONSOON with this tree. */
    Queue<Bit> bits = {};

    EXPECT_EQUAL(decodeText(bits, tree), "");

    deleteTree(tree);
}


ADD_TEST("Provided Test: encodeText returns an empty queue when given an empty input.") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N
     *                / \
     *               M   S
     */
    HuffmanNode* tree = new HuffmanNode {
           '*',
           new HuffmanNode { 'O', nullptr, nullptr },
               new HuffmanNode { '*',
               new HuffmanNode{ '*',
               new HuffmanNode { 'M', nullptr, nullptr },
               new HuffmanNode { 'S', nullptr, nullptr }
            },
            new HuffmanNode{ 'N', nullptr, nullptr }
        }
    };

    /* What you get if you encode MONSOON with this tree. */
    Queue<Bit> expected = {};

    EXPECT_EQUAL(encodeText("", tree), expected);

    deleteTree(tree);
}


ADD_TEST("Provided Test: encode can handle a tree with two leaf nodes of equal weight.") {

    HuffmanNode* tree = huffmanTreeFor("AB");

    Queue<Bit>  bits;
    Queue<char> leaves;

    encodeTree(tree, bits, leaves);

    Queue<Bit>  expectedBits   = { 1, 0, 0 };
    Queue<char> expectedLeaves = { 'A', 'B' };

    EXPECT_EQUAL(bits,   expectedBits);
    EXPECT_EQUAL(leaves, expectedLeaves);

    deleteTree(tree);
}


bool isEncodingTree(HuffmanNode* tree);
bool areEqual(HuffmanNode* lhs, HuffmanNode* rhs);
ADD_TEST("Provided Test: Can decode a tree with a single leaf node.") {

    Queue<Bit>  bits   = { 0 };
    Queue<char> leaves = { 'A' };

    HuffmanNode* tree = decodeTree(bits, leaves);
    EXPECT(isEncodingTree(tree));

    /* Confirm this is the right tree. */
    HuffmanNode* expected = new HuffmanNode { 'A', nullptr, nullptr };

    EXPECT(areEqual(tree, expected));

    deleteTree(tree);
    deleteTree(expected);
}







/* * * * * Provided Tests Below This Point * * * * */
#include <limits>

/* Utility function to test if a purported Huffman tree is indeed a Huffman tree.
 * Specifically, this checks that each internal node has either zero or two
 * children. There are other ways you could produce an invalid Huffman tree - for
 * example, by having uninitialized pointers or by linking in a cycle - but we
 * don't test for that here.
 */
bool isEncodingTree(HuffmanNode* tree) {
    /* The empty tree is not a Huffman tree. */
    if (tree == nullptr) return false;

    /* If we have one missing child, we should have two missing children. */
    if ((tree->zero == nullptr) != (tree->one == nullptr)) return false;

    /* If we have children at all, they need to be Huffman trees. */
    return tree->zero == nullptr || (isEncodingTree(tree->zero) && isEncodingTree(tree->one));
}

/* Utility function to test if two trees are equal. This is adapted from Section
 * Handout 8 and particularized to Huffman trees.
 */
bool areEqual(HuffmanNode* lhs, HuffmanNode* rhs) {
    /* Base case: If either is a leaf, both should be. */
    bool lhsLeaf = lhs->zero == nullptr && lhs->one == nullptr;
    bool rhsLeaf = rhs->zero == nullptr && rhs->one == nullptr;
    if (lhsLeaf || rhsLeaf) {
        return lhs->ch == rhs->ch && lhsLeaf == rhsLeaf;
    }

    /* Otherwise, they're both internal nodes. Check that they match. */
    return areEqual(lhs->zero, rhs->zero) && areEqual(lhs->one, rhs->one);
}

/* Utility function to return a string of all possible characters. */
string pangrammaticString() {
    string result;

    char ch = numeric_limits<char>::min();
    result += ch;
    do {
        ch++;
        result += ch;
    } while (ch != numeric_limits<char>::max());

    return result;
}

/* Utility function that makes an inefficient but still valid encoding tree
 * for the given characters.
 */
HuffmanNode* strandTreeFor(const string& text, size_t index = 0) {
    if (index == text.size()) error("No characters provided to strandTreeFor.");

    /* We always get a leaf node. */
    HuffmanNode* leaf = new HuffmanNode {
        text[index], nullptr, nullptr
    };

    /* Last character? If so, that's all. */
    if (index + 1 == text.size()) return leaf;

    /* Otherwise, build a larger tree. */
    else return new HuffmanNode {
        ' ', leaf, strandTreeFor(text, index + 1)
    };
}

ADD_TEST("Provided Test: huffmanTreeFor reports errors on invalid inputs.") {
    EXPECT_ERROR(huffmanTreeFor(""));    // No characters
    EXPECT_ERROR(huffmanTreeFor("a"));   // One character
    EXPECT_ERROR(huffmanTreeFor("aaa")); // One character
}

ADD_TEST("Provided Test: huffmanTreeFor builds tree for two characters.") {
    HuffmanNode* reference = new HuffmanNode {
        ' ', new HuffmanNode {'a', nullptr, nullptr}, new HuffmanNode {'b', nullptr, nullptr}
    };

    HuffmanNode* tree = huffmanTreeFor("aaabbbb");
    EXPECT(isEncodingTree(tree));
    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

ADD_TEST("Provided Test: huffmanTreeFor works on the full range of characters.") {
    /* Get a string of all possible characters, then pair them off and see what we find. */
    string allChars = pangrammaticString();
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string toEncode;
        toEncode += allChars[i];
        toEncode += allChars[i + 1];
        toEncode += allChars[i + 1];

        HuffmanNode* reference = new HuffmanNode {
            ' ',
            new HuffmanNode {allChars[i], nullptr, nullptr},
            new HuffmanNode {allChars[i + 1], nullptr, nullptr}
        };

        HuffmanNode* tree = huffmanTreeFor(toEncode);
        EXPECT(isEncodingTree(tree));
        EXPECT(areEqual(tree, reference));

        deleteTree(reference);
        deleteTree(tree);
    }
}

ADD_TEST("Provided Test: huffmanTreeFor uses cumulative weights.") {
    /* This tree:
     *                 *
     *                / \
     *               *   D
     *              / \
     *             C   *
     *                / \
     *               A   B
     */
    HuffmanNode* reference = new HuffmanNode {
        '*',
            new HuffmanNode { '!',
                new HuffmanNode { 'C', nullptr, nullptr },
                new HuffmanNode { '?',
                    new HuffmanNode { 'A', nullptr, nullptr },
                    new HuffmanNode { 'B', nullptr, nullptr }
                }
            },
            new HuffmanNode { 'D', nullptr, nullptr }
    };

    /* Ax2, Bx3, Cx4, Dx10 */
    HuffmanNode* tree = huffmanTreeFor("AABBBCCCCDDDDDDDDDD");
    EXPECT(isEncodingTree(tree));
    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

ADD_TEST("Provided Test: decodeText works on small sample.") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N
     *                / \
     *               M   S
     */
    HuffmanNode* tree = new HuffmanNode {
        '*',
            new HuffmanNode { 'O', nullptr, nullptr },
            new HuffmanNode { '*',
                new HuffmanNode{ '*',
                    new HuffmanNode { 'M', nullptr, nullptr },
                    new HuffmanNode { 'S', nullptr, nullptr }
                },
                new HuffmanNode{ 'N', nullptr, nullptr }
            }
    };

    /* What you get if you encode MONSOON with this tree. */
    Queue<Bit> bits = { 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1 };

    EXPECT_EQUAL(decodeText(bits, tree), "MONSOON");

    deleteTree(tree);
}

ADD_TEST("Provided Test: Can decode all char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try decoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string expected;
        expected += allChars[i];
        expected += allChars[i + 1];
        expected += allChars[i + 1];

        HuffmanNode* tree = new HuffmanNode {
            ' ',
            new HuffmanNode {allChars[i], nullptr, nullptr},
            new HuffmanNode {allChars[i + 1], nullptr, nullptr}
        };

        /* Decode the bitstream 011, which should map back to the expected
         * string.
         */
        Queue<Bit> bits = { 0, 1, 1 };
        EXPECT_EQUAL(decodeText(bits, tree), expected);

        deleteTree(tree);
    }
}

ADD_TEST("Provided Test: encodeText works on small sample.") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N
     *                / \
     *               M   S
     */
    HuffmanNode* tree = new HuffmanNode {
           '*',
           new HuffmanNode { 'O', nullptr, nullptr },
               new HuffmanNode { '*',
               new HuffmanNode{ '*',
               new HuffmanNode { 'M', nullptr, nullptr },
               new HuffmanNode { 'S', nullptr, nullptr }
            },
            new HuffmanNode{ 'N', nullptr, nullptr }
        }
    };

    /* What you get if you encode MONSOON with this tree. */
    Queue<Bit> expected = { 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1 };

    EXPECT_EQUAL(encodeText("MONSOON", tree), expected);

    deleteTree(tree);
}

ADD_TEST("Provided Test: Can encode all char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string toEncode;
        toEncode += allChars[i];
        toEncode += allChars[i + 1];
        toEncode += allChars[i + 1];

        HuffmanNode* tree = new HuffmanNode {
                ' ',
                new HuffmanNode {allChars[i], nullptr, nullptr},
                new HuffmanNode {allChars[i + 1], nullptr, nullptr}
        };

        /* See what bits we get back. We should get 011, since the first
         * character has code 0 and the second has code 1.
         */
        Queue<Bit> bits = encodeText(toEncode, tree);
        Queue<Bit> expected = { 0, 1, 1 };

        EXPECT_EQUAL(bits, expected);

        deleteTree(tree);
    }
}

ADD_TEST("Provided Test: decodeText undoes encodeText on range of sample strings.") {
    Vector<string> testCases = {
        "THAT THAT IS IS THAT THAT IS NOT IS NOT IS THAT IT IT IS",
        "AABAAABBABAAABAAAA",
        ":-) :-D XD <(^_^)>",
        pangrammaticString(),
    };

    for (string test: testCases) {
        /* Use a silly encoding scheme, but one that works regardless of what
         * characters are provided.
         */
        HuffmanNode* tree = strandTreeFor(test);
        EXPECT(isEncodingTree(tree));

        Queue<Bit> bits = encodeText(test, tree);
        string result = decodeText(bits, tree);

        /* Check if the strings are equal. Don't use EXPECT_EQUAL so that we
         * don't blow up the console if there's a mismatch.
         */
        EXPECT_EQUAL(test.size(), result.size());
        EXPECT(test == result);

        deleteTree(tree);
    }
}

ADD_TEST("Provided Test: Can encode an example tree.") {
    /* Build an encoding tree for "ABBCCCC." It should look like this:
     *
     *                 *
     *                / \
     *               *   C
     *              / \
     *             A   B
     *
     * This will compress down to
     *
     *        11000
     *        ABC
     */
    HuffmanNode* tree = huffmanTreeFor("ABBCCCC");

    Queue<Bit>  bits;
    Queue<char> leaves;

    encodeTree(tree, bits, leaves);

    Queue<Bit>  expectedBits   = { 1, 1, 0, 0, 0 };
    Queue<char> expectedLeaves = { 'A', 'B', 'C' };

    EXPECT_EQUAL(bits,   expectedBits);
    EXPECT_EQUAL(leaves, expectedLeaves);

    deleteTree(tree);
}

ADD_TEST("Provided Test: Can decode an example tree.") {
    /* This encodes this tree:
     *
     *                 *
     *                / \
     *               *   C
     *              / \
     *             A   B
     */
    Queue<Bit>  bits   = { 1, 1, 0, 0, 0 };
    Queue<char> leaves = { 'A', 'B', 'C' };

    HuffmanNode* tree = decodeTree(bits, leaves);
    EXPECT(isEncodingTree(tree));

    /* Confirm this is the right tree. */
    HuffmanNode* expected = new HuffmanNode {
        '*',
            new HuffmanNode {
                '*',
                new HuffmanNode { 'A', nullptr, nullptr },
                new HuffmanNode { 'B', nullptr, nullptr },
            },
            new HuffmanNode { 'C', nullptr, nullptr }
    };

    EXPECT(areEqual(tree, expected));

    deleteTree(tree);
    deleteTree(expected);
}

ADD_TEST("Provided Test: Can encode trees using all possible char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        HuffmanNode* tree = new HuffmanNode {
            ' ',
            new HuffmanNode {allChars[i], nullptr, nullptr},
            new HuffmanNode {allChars[i + 1], nullptr, nullptr}
        };

        /* See what we get back. It should be the bitstring 100 (root with
         * two children) and the two leaves, in order.
         */
        Queue<Bit>  treeBits;
        Queue<char> treeLeaves;

        Queue<Bit>  expectedBits = { 1, 0, 0 };
        Queue<char> expectedLeaves = { allChars[i], allChars[i + 1] };

        encodeTree(tree, treeBits, treeLeaves);
        EXPECT_EQUAL(treeBits, expectedBits);
        EXPECT_EQUAL(treeLeaves, expectedLeaves);

        deleteTree(tree);
    }
}

ADD_TEST("Provided Test: Can decode trees using all possible char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        HuffmanNode* expected = new HuffmanNode {
            ' ',
            new HuffmanNode {allChars[i], nullptr, nullptr},
            new HuffmanNode {allChars[i + 1], nullptr, nullptr}
        };
        Queue<Bit>  treeBits   = { 1, 0, 0 };
        Queue<char> treeLeaves = { allChars[i], allChars[i + 1] };

        HuffmanNode* tree = decodeTree(treeBits, treeLeaves);
        EXPECT(isEncodingTree(tree));
        EXPECT(areEqual(tree, expected));

        deleteTree(tree);
        deleteTree(expected);
    }
}

ADD_TEST("Provided Test: decodeTree undoes encodeTree on sample strings.") {
    /* Make a Huffman tree for the string of all characters. */
    HuffmanNode* sourceTree = huffmanTreeFor(pangrammaticString());
    EXPECT(isEncodingTree(sourceTree));

    /* Encode, then decode it. */
    Queue<Bit>  bits;
    Queue<char> leaves;
    encodeTree(sourceTree, bits, leaves);

    HuffmanNode* resultTree = decodeTree(bits, leaves);
    EXPECT(isEncodingTree(resultTree));
    EXPECT(areEqual(sourceTree, resultTree));

    deleteTree(sourceTree);
    deleteTree(resultTree);
}

ADD_TEST("Provided Test: Can decompress a small sample file.") {
    HuffmanResult file = {
        { 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0 },
        { 'u', 'k', 'p', 'n', 'a', 'm', 'h' },
        { 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1,
          0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0,
          0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0 }
    };

    EXPECT_EQUAL(decompress(file), "humuhumunukunukuapuaa");
}

ADD_TEST("Provided Test: Compress reports errors on bad inputs.") {
    EXPECT_ERROR(compress(""));
    EXPECT_ERROR(compress("A"));
    EXPECT_ERROR(compress("AAAA"));
}

ADD_TEST("Provided Test: Can compress a small sample file.") {
    HuffmanResult file = compress("ABANANAABANDANA");
    Queue<Bit>  treeBits    = { 1, 1, 1, 0, 0, 0, 0 };
    Queue<char> treeChars   = { 'D', 'B', 'N', 'A' };
    Queue<Bit>  messageBits = { 1, 0, 0, 1, 1, 0, 1, 1, 0,
                                1, 1, 1, 0, 0, 1, 1, 0, 1,
                                0, 0, 0, 1, 0, 1, 1 };

    EXPECT_EQUAL(file.treeBits, treeBits);
    EXPECT_EQUAL(file.treeLeaves, treeChars);
    EXPECT_EQUAL(file.messageBits, messageBits);
}

ADD_TEST("Provided Test: Compress undoes decompress on a range of strings.") {
    Vector<string> testCases = {
        "THAT THAT IS IS THAT THAT IS NOT IS NOT IS THAT IT IT IS",
        "AABAAABBABAAABAAAA",
        ":-) :-D XD <(^_^)>",
        pangrammaticString(),
    };

    for (string test: testCases) {
        HuffmanResult file = compress(test);
        string result = decompress(file);

        EXPECT_EQUAL(result.size(), test.size());

        /* Don't clobber the output with a huge string if there's a mismatch. */
        EXPECT(test == result);
    }
}
