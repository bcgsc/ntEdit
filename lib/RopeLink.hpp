/*
 *
 * RopeLink Data Structure
 * Author: Jessica Zhang
 * Genome Sciences Centre, 
 * British Columbia Cancer Agency 
 */

#ifndef ROPELINK_H_
#define ROPELINK_H_
#endif

#include "iostream" 
#include <cstdlib>
#include <string> 

using namespace std; 

class RopeLink {
	public: 
		RopeLink *left, *right; // [left]---[this]---[right]
		int node_type; // 0 = position node; 1 = character node; 2 = string node;
		size_t s_pos; // 0-based
		size_t e_pos; // 0-based
		unsigned char c;
		char *str;  

	// creates a RopeLink that holds a string
//	void createStringRopeLink(RopeLink *&node, RopeLink *&left, string& to_store) {
//		RopeLink *tmp = new RopeLink(); 
//		tmp->left = left; 
//		tmp->right = NULL; 
//		tmp->str = to_store.c_str(); 
//
//		node = tmp; 
//	}

	RopeLink() {
	}

	void createCharacterRopeLink(RopeLink *&node, RopeLink *&left, unsigned char c) {
		RopeLink *tmp = new RopeLink(); 
		tmp->left = left; 
		tmp->right = NULL; 
		tmp->str = NULL; 
		tmp->c = c; 
		
		node=tmp; 
	}
	
	// creates a RopeLink that holds positions
	void createPositionRopeLink(RopeLink *&node, RopeLink *&left, size_t s_pos, size_t e_pos) {
		RopeLink *tmp = new RopeLink(); 
		tmp->left = left; 
		tmp->right = NULL; 
		tmp->str = NULL; 
		tmp->s_pos = s_pos; 
		tmp->e_pos = e_pos; 

		node = tmp; 
	}

	// make an insertion with a string before a specified position within a position rope node
	// 	example: if I make an insertion at position 5: 
	// 			ABCDEFG then insert aaa at position 5
	// 			result: ABCDEaaaFG
//	void makeStringInsertion(RopeLink *&node, size_t insert_pos, string& to_insert) {
//		RopeLink *tmp_right_node = node->right;
//		size_t tmp_epos = node->e_pos;
//		node->e_pos = insert_pos-1; // assumes that insert_pos > 0
//		RopeLink *insertion_node = NULL; 
//		createStringRopeLink(insertion_node, node, to_insert); 
//		node->right = insertion_node; 
//		RopeLink *right_node = NULL; 
//		createPositionRopeLink(right_node, insertion_node, insert_pos, tmp_epos); 
//		insertion_node->right = right_node;
//		right_node->right = tmp_right_node; 
//	}

	// make a character insertion into a string right before a specified position in a position rope node
	// 	example: insert 'a' into ABCDEFG at position 3
	// 	result: ABC a DEFG
	void makeCharacterInsertionOnPosition(RopeLink *&node, size_t insert_pos, unsigned char c) {
		RopeLink *tmp_right_node = node->right;
		size_t tmp_epos = node->e_pos; 
		node->e_pos = insert_pos-1; 
		RopeLink *insertion_node = NULL; 
		createCharacterRopeLink(insertion_node, node, c); 
		node->right = insertion_node; 
		RopeLink *right_node = NULL; 
		createPositionRopeLink(right_node, insertion_node, insert_pos, tmp_epos); 
		insertion_node->right = right_node;
		right_node->right = tmp_right_node;
	}

	// make a character insertion after a node that is also a character or string RopeLink
	void makeCharacterInsertionOnCharacter(RopeLink *&node, unsigned char c) {
		RopeLink *tmp_right_node = node->right; 
		RopeLink *insertion_node = NULL; 
		createCharacterRopeLink(insertion_node, node, c); 
		node->right = insertion_node; 
		insertion_node->right = tmp_right_node; 
	}

	void makeCharacterInsertion(RopeLink *&node, int insert_pos, unsigned char c) {
		if (node->node_type == 0) 
			makeCharacterInsertionOnPosition(node, insert_pos, c); 
		else if (node->node_type == 1) 
			makeCharacterInsertionOnCharacter(node, c); 
	}


	// make a deletion with a string including the specified position of a position rope node
	// 	example: if I make a deletion at position 3 of 2 characters: 
	// 			ABCDEFG delete 2 at position 3: 
	// 			result: ABC	FG
	void makeDeletion(RopeLink *&node, size_t delete_pos, size_t num_del) {
		RopeLink *tmp_right_node = node->right; 
		size_t tmp_epos = node->e_pos; 
		node->e_pos = delete_pos-1; 
		RopeLink *right_node = NULL; 
		createPositionRopeLink(right_node, node, delete_pos+num_del, tmp_epos); 
		node->right = right_node;
		right_node->right = tmp_right_node; 
	}

	// Remove a specified node
	void removeRopeLink(RopeLink *& to_remove) {
		RopeLink * prev_node = to_remove->left;
		prev_node->right = to_remove->right;
		RopeLink * next_node = to_remove->right;
		next_node->left = prev_node;

		to_remove->left = NULL; 
		to_remove->right = NULL;
		delete to_remove; 
	}

	/** destructor */
	~RopeLink() {
		if (str!=NULL)
			delete [] str; 
	}
		

};
