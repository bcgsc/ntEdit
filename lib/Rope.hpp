/*
 *
 * Rope Data Structure
 * Author: Jessica Zhang
 * Genome Sciences Centre, 
 * British Columbia Cancer Agency 
 */

#ifndef ROPE_H_
#define ROPE_H_

#include "iostream" 
#include <cstdlib>
#include <string> 

using namespace std; 

class Rope {
	public: 
		Rope *left, *right, *parent; 
		char *str; 
	}; 

	// creates a Rope node that holds a string
	void createRope(Rope *&node, Rope *&parent, string& to_store) {
		Rope *tmp = new Rope(); 
		tmp->left = tmp->right = NULL; 
		tmp->parent = parent; 

		node = tmp; 
		// TODO: check that it is okay that we set node before AND that we stored a reference to a string to_store
		tmp->str = to_store.c_str(); 
	}

	// links two strings to a parent
	Rope * concatenateTwoRopes(Rope *&word1, Rope*&word2) {
		Rope *tmp = new Rope(); 
		tmp->left = word1; 
		tmp->right = word2; 
		word1->parent = tmp; 
		word2->parent = tmp; 

		tmp->str = NULL; 

		return tmp; 
	}

	// insert a string to_insert after a word
	void insertRope(Rope *&word, Rope*& to_insert) {
		Rope * tmp = new Rope(); 
		tmp->parent = word->parent; 
		tmp->left = word; 
		tmp->right = to_insert; 
		word->parent = tmp; 
		to_insert->parent = tmp; 

		tmp->str = NULL; 

	}

	// read the full string using inorder traversal 
//	string readFullRope(Rope *& root, unsigned est_size) {
//		string fullSequence; 
//		fullSequence.reserve(est_size); 
//		stringstream ss(fullSequence); 
//
//		stack<Rope *> todo;
//		Rope * curr = root; 
//		while (curr != NULL || !todo.empty()) {
//			while (curr != NULL) {
//				todo.push(curr); 
//				curr = curr->left; 
//			}
//			curr = todo.top(); 
//			todo.pop(); 
//
//			if (curr->str != NULL) {
//				ss << curr->str; 
//			}
//
//			curr = curr->right; 
//		}
//		return ss; 
//	}
};
