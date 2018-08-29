#include <iostream>
#include <thread>
#include <mutex>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
using namespace std;

mutex all;

class Linked_list {

	struct node {
		// Initialise a mutex to lock the node
		mutex m;
		// Id string for the node
		string Identifier_string;
		//Each node constrained by the nodes before and after it
		node* next;
		node* prev;
	};
	// Define the beginnig and end of double linked list
	node* head;
	node* tail;

	public:
		// Value of the current length of the list
		int Current_length;
		// Generate a random string for the node (consisting 2-8 characters)
		string Generate_string();
		// Create a new node
		void Create_node();
		// Remove a random node from the list
		void Remove_Random_Node();
		// Traverse the list and concatenate all strings of the nodes
		void Traverse_list();
};

int main(int argc, char** argv) {
	//Adding 100 nodes to the doubly linked list
	// Initialise random seed for number generator
	srand(time(0));
	// Define the number of nodes in the list to be created
	int List_size = 100;
	// Create a new instance of the Linked list
	Linked_list list = Linked_list();
	// Keep appending new nodes to the list until it reaches required size
	for (int i=0; i<List_size; i++) {
		list.Create_node();
	}

	// First thread to traverse the list and concatenate node strings until no nodes left
	// The thread uses a function from the main thread by reference
	thread Traversing_thread([&]() {
		while (list.Current_length != 0) {
			// Lock other threads so concatenated string can be printed without being interrupted
			all.lock();
			cout << "The Concatenated node string is:"<<endl;
			list.Traverse_list();
			all.unlock();
		}
	});

	// Second thread to keep deleting random nodes until the list is empty
	thread Removing_thread([&]() {
		while (list.Current_length > 0) {
			list.Remove_Random_Node();
			// Delay 1s before removing next node
			this_thread::sleep_for(chrono::milliseconds(1000)); 
		}
	});

	// Join both threads to ensure they finish executing and run in the background
	Traversing_thread.join();
	Removing_thread.join();

	return 0;
}


string Linked_list :: Generate_string() {
	//Dictionary of letters to be used
	static const char Letters[] =
	"abcdefghijklmnopqrstuvwxyz";

	int Number_of_letters = sizeof(Letters) - 1;
	string Random_string;

	//generate random string length in the range 2-8 inclusive
	int Random_length = rand() % 7 + 2;
	for(int i = 0; i<Random_length; i++) {
		//append string with random letter
	 Random_string += Letters[rand() % Number_of_letters];
	}
	return Random_string;
}


void Linked_list :: Create_node() {
	// Create a new node n inside our List
	node* n = new node;
	// Assign the new node a random string
	n->Identifier_string = Generate_string();

	// If the list is empty, initialise the first node
	if (Current_length == 0){
		// First node will not have any nodes in front or behind it
		n->prev = NULL;
		n->next = NULL;
		// First node will be the head and tail of the list at the same time
		head = n;
		tail = n;
	}

	// Append a new node if list isnt empty
	else {
		// Lock the node's tail before defining its parents and children
		lock_guard<mutex> Tail_Lock(tail->m);
		// Append the new node at the end of the list
		// Set new node's parent as the previous end of the list
		n->prev = tail;
		// The previous end of the list now has the new node as its child
		tail->next = n;
		// Set the new node as the new end of the list
		tail = n;
	}
	// Increment the list's length after adding the new node
	Current_length = Current_length + 1;
}


void Linked_list :: Remove_Random_Node() {
	//Only attempt to remove nodes when the list is not empty
	while (Current_length != 0) {
		// Lock the head node before setting a pointer to it
		unique_lock<mutex> Unique_Lock(head->m);
		node* Selected_node = head;
		Unique_Lock.unlock();

		// Locate a randomly positioned node
		int Random_position = rand() % Current_length;
		// Increment the head node until the random position
		for(int i = 0; i<Random_position; i++) {
			// Lock both the current node and the one ahead of it
			lock_guard<mutex> Node_Lock(Selected_node->m);
			lock_guard<mutex> Next_Node_Lock(Selected_node->next->m);
			Selected_node = Selected_node->next;
		}
		// Lock everything before printing
		all.lock();
		cout <<"Removing node at position = " << Random_position+1<<", List Length = " << Current_length<<endl;
		all.unlock();

		// Handling different node deletion cases
		// Locking order should be from left to right to avoid deadlock
		// If only one node left in the list, remove both head and tail pointers
		if (Current_length == 1) {
			lock_guard<mutex> Node_Lock(Selected_node->m);
			head = NULL;
			tail = NULL;
			Selected_node = NULL;
		}
		// If the first node (head) is to be deleted
		else if (Selected_node == head) {
			// Lock current node first, then the next one
			lock_guard<mutex> Node_Lock(Selected_node->m);
			lock_guard<mutex> Next_Node_Lock(Selected_node->next->m);
			// The new head is the second node in the list
			head = head->next;
			// Head node no has no parent
			head->prev = NULL;
		}
		// If the last node (tail) is to be deleted
		else if (Selected_node == tail) {
			// Lock the previous node first, then the current one
			lock_guard<mutex> Previous_Node_Lock(Selected_node->prev->m);
			lock_guard<mutex> Node_Lock(Selected_node->m);
			// The new tail is the second last node
			tail = Selected_node->prev;
			// Tail node has no children
			tail->next = NULL;
		}
		// Otherwise, if the node is in the middle of the list
		else {
			// Lock first the previous node, then the current and then the next
			lock_guard<mutex> Previous_Node_Lock(Selected_node->prev->m);
			lock_guard<mutex> Node_Lock(Selected_node->m);
			lock_guard<mutex> Next_Node_Lock(Selected_node->next->m);
			// Node behind the selected node must point to the node in front of the selected node
			Selected_node->prev->next = Selected_node->next;
			// Node in front of the selected node must point to the node before the selected node
			Selected_node->next->prev = Selected_node->prev;
		}
		// Delete the selected node from memory after handling it
		free(Selected_node);
		// Decrement list length after node deletion
		Current_length = Current_length - 1;
	}
}


void Linked_list :: Traverse_list () {
	// Ensure that the list isn't empty
	if (Current_length != 0) {
		// Store the concatenated string
		string Concatenated_string;
		// Lock the head node before assigning a pointer to it
		unique_lock<mutex> Unique_Lock(head->m);
		// Set the current node pointer to head (beginning of the list)
		node* Selected_node = head;
		Unique_Lock.unlock();
		// Loop through the list, concatenating all node strings, except for the last one
		// Attempting to traverse a single node list results in segmentation faults
		for (int i = 0; i<Current_length-1; i++){
			// Lock both the current node and the one ahead of it
			lock_guard<mutex> Node_Lock(Selected_node->m);
			lock_guard<mutex> Next_Node_Lock(Selected_node->next->m);
			// Append the current node string 
			Concatenated_string += Selected_node->Identifier_string;
			// Go to the next node
			Selected_node = Selected_node->next;
		}
		//Append string of the last node
		Concatenated_string += Selected_node->Identifier_string;
		cout<<Concatenated_string<<endl;
	}
}