#ifndef MINIMAL_LIST_H
#define MINIMAL_LIST_H

// Minimal list to store structs that otherwise need annoying deleted constructors
template <typename T>
struct minimal_list {
	struct node {
		T value;
		node* next = 0;
	};

	minimal_list() = default;
	~minimal_list() {
		if (!head)
			return;
		node* curr = head;
		node* next = head->next;
		while (next) {
			delete curr;
			curr = next;
			next = next->next;
		}
	}

	void Push(const T& val) {

		if (!head) {
			head = new node;
			head->value = val;
		}
		else {
			node* end = head;
			while (end->next) {
				end = end->next;
			}
			end->next = new node;
			end->next->value = val;
		}
		++size;
	}

	// Get entry at index i (assumes that i is within bounds)
	node* Entry(int i) {
		node* e = head;
		for (int it = 0; it < i; it++)
			e = head->next;

		return e;
	}

	node* End() {
		return Entry(size - 1);
	}

	void Pop() {

		if (!head)
			return;

		node* curr = head;
		node* next = head->next;

		// Only head
		if (!next) {
			delete head;
			head = 0;
			--size;
			return;
		}

		while (next) {

			// next is the end of the list
			// so delete next and set curr->next to zero
			if (!next->next)
				break;

			curr = next;
			next = next->next;
		}

		delete next;
		curr->next = 0;
		--size;
	}


	int Size() const { return size; }


	node* head = 0;
	int size = 0; // Default list size is 0
};


#endif