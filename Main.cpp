/*
* Main.cpp
*
*  Created on: Sep 12, 2017
*      Author: julio
**/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <functional>

using namespace std;

string fileNames[] = {
	"Almost100.txt","Almost1000.txt",
	"Inverse100.txt","Inverse1000.txt",
	"Random100.txt","Random1000.txt"
};


//Create 6 arrays for the 6 files and initialize them to 0
int Almost100[100] = { 0 },
Almost1000[1000] = { 0 },
Inverse100[100] = { 0 },
Inverse1000[1000] = { 0 },
Random100[100] = { 0 },
Random1000[1000] = { 0 };

//Create a pointer that point  to the 6 arrays. 
int *bigArray[6] = {
	Almost100,Almost1000,
	Inverse100,Inverse1000,
	Random100,Random1000
};

int comparisons = 0;
int swaps = 0;


void selectionSort(int array[], int size);
void insertionSort(int array[], int size);
void mergeSort(int array[], int l, int r);
void merge(int array[], int l, int m, int r);
void heapSort(int array[], int size);
void heapify(int arr[], int n, int i);
void quickSort(int array[], int low, int high);
int partition(int array[], int low, int high);
void countingSort(int array[], int size, int range);
void swap(int&, int&);

void fillArray(int**,int size);                             //Fill an array of arrays.
void getInputFrom(ifstream& file, int array[], int size);	//Get input from Files.
void printData(string fileName, int swaps, int comparisons); //Print the number of comparisons and swaps
void printHeader(string sortName);                           //Create a tittle
void analyzeAlgorithm(void(*function)(int*, int));           //Analyze an algorithm by taking in a function
void analyzeAlgorithm(void(*function)(int*, int, int));     //Analyze an algorithm by taking in a function.
int main()
{
	//Names of the sorting algorithms
	string sortNames[] = {
		"Selection","Insertion",
		"Merge","heap","quick",
		"Counting"
	};

	//selectionSort
	fillArray(bigArray, 6);
	printHeader(sortNames[0]);
	analyzeAlgorithm(selectionSort);

	cout << endl;
	//InsertionSort
	fillArray(bigArray, 6);
	printHeader(sortNames[1]);
	analyzeAlgorithm(insertionSort);

	cout << endl;
	//Mergesort
	fillArray(bigArray, 6);
	printHeader(sortNames[2]);
	analyzeAlgorithm(mergeSort);

	cout << endl;
	//HeapSort
	fillArray(bigArray, 6);
	printHeader(sortNames[3]);
	analyzeAlgorithm(heapSort);

	cout << endl;
	//QuickSort
	fillArray(bigArray, 6);
	printHeader(sortNames[4]);
	analyzeAlgorithm(quickSort);

	
	cout << endl;
	//QuickSort
	fillArray(bigArray, 6);
	printHeader(sortNames[5]);
	bool flag = false;
	for (int i = 0; i < 6; i++) {
		int n = flag ? 1000 : 100;
		countingSort(bigArray[i], n, n + 1);
		printData(fileNames[i], swaps, comparisons);
		comparisons = 0;
		swaps = 0;
		flag = !flag;
	}

	//Prevent visual studio to close output
	int wait = 0;
	cin >> wait;
	return 0;
}

void analyzeAlgorithm(void(*function)(int*, int, int)){
	bool flag = false;
	for (int i = 0; i < 6; i++) {
		int n = flag ? 1000 : 100;
		(*function)(bigArray[i], 0,n);
		printData(fileNames[i], swaps, comparisons);
		comparisons = 0;
		swaps = 0;
		flag = !flag;
	}
}
void analyzeAlgorithm(void(*function)(int*, int)) {

	bool flag = false;
	for (int i = 0; i < 6; i++){
		int n = flag ? 1000 : 100;
		(*function)(bigArray[i],n);
		printData(fileNames[i], swaps, comparisons);
		comparisons = 0;
		swaps = 0;
		flag = !flag;
	}
}

/*************************************************************
|The printHeader function recives the name of the algorithm 
|being analize and create a header format to display.
**************************************************************/
void printHeader(string sortName){
	cout << setw(20) << sortName << setw(20)
		<< "Swaps" << setw(20) << "Comparisons" << endl;

	cout << setw(60) << "________________________________________________________" << endl;
}
/****************************************************************
|The printData function receives the name of the file analized,
|the number of swaps, and comparisons and display it on the
|screen.
****************************************************************/
void printData(string fileName, int swaps, int comparisons) {
	cout  << setw(20) << fileName << setw(20) 
		<< swaps  << setw(20) << comparisons << endl;
}

/************************************************
|The fillArray function  fill an array of arrays
|by reading data from 6 files and putting it in.
************************************************/
void fillArray(int** array,int size){

	ifstream file[] = {
		ifstream("Almost100.txt"),ifstream("Almost1000.txt"),
		ifstream("Inverse100.txt"),ifstream("Inverse1000.txt"),
		ifstream("Random100.txt"),ifstream("Random1000.txt")
	};

	bool flag = false;
	for (int i = 0; i < size; i++)
	{
		if (!flag){
			getInputFrom(file[i], array[i], 100);
			flag = true;
		}
		else {
			getInputFrom(file[i], array[i], 1000);
			flag = false;
		}
	}
}

/**********************************************
|The getInputFrom function receives a file,
|read the data,store it in an array of integers
|and return it.
***********************************************/
void getInputFrom(ifstream& inputFile, int array[], int size)
{
	for (int i = 0; i < size; i++)
		inputFile >> array[i];
	inputFile.close();
}

/**********************************************
|The swap function, receives two integers
|passed by reference and swap them.
***********************************************/
void swap(int&value1, int&value2) {
	int temp = value1;
	value1 = value2;
	value2 = temp;
}

/************************************************
|The selectionSort function receives an array of
|integers and its size and sort the array by
|finding the smallest element repeatedly and
|placing it in decending order.
************************************************/
void selectionSort(int array[], int size) {
	for (int i = 0; i < size - 1; i++) {
		int min = i;
		for (int j = i + 1; j < size; j++) {
			if (array[j]<array[min])
				min = j;
			comparisons++;
		}
		swaps++;
		swap(array[min], array[i]);
	}
}

/********************************************************
|The insertionSort function receives an array of
|integers and its size and sort the array by making
|the given array into two sublist the sorted and unsorted.
|The sorted starting with its first element keep appending
|element into the sorted list until the unsorted is empty.
***********************************************************/
void insertionSort(int array[], int size) {
	int i, key, j;
	for (i = 1; i < size; i++){
		key = array[i];
		j = i - 1;

		while (j >= 0 && array[j] > key){
			array[j + 1] = array[j];
			j = j - 1;
			swaps++;
			comparisons++;
		}
		array[j + 1] = key;
		comparisons++;
	}
}

/**************************************************************
|The mergeSort function receives an array of integers and
|the beginning and end of the array as the left and right
|sides. The array it sorted by dividing the array into sublists
|until each sublists contain one element, and merge them
|together back in order by calling the merge funcition.
**************************************************************/
void mergeSort(int array[], int l, int r) {
	
	if (l < r) {
		comparisons++;
		int middle = (l + r) / 2;

		mergeSort(array, l, middle);
		mergeSort(array, middle + 1, r);

		merge(array, l, middle, r);
	}
}

/*****************************************************
|The merge function receives an array with two sorted
|sublists and merge them together.
|The l parameter is the beginning and m is the last
|element of the first sublist. m+1 is the biginning and
|r is the last element in the second sublist.
******************************************************/
void merge(int array[], int l, int m, int r) {
	int i, j, k;

	//Get the length of  left sublist
	int n1 = m - l + 1;

	//Get the length of the right sublist
	int n2 = r - m;

	/*create temp arrays*/
	int* L = new int[n1];
	int* R = new int[n2];

	//Copy data to temp arrays L[] and R[]
	for (i = 0; i < n1; i++)
		L[i] = array[l + i];
	for (j = 0; j < n2; j++)
		R[j] = array[(m + 1) + j];

	i = 0; 
	j = 0; 
	k = l;
	while (i < n1 && j < n2) {
		if (L[i] <= R[j])
			array[k] = L[i++];
		else
			array[k] = R[j++];
		k++;
		swaps++;
		comparisons++;
	}

	/*Copy the remaining elements of L[], if there are any*/
	while (i < n1){
		array[k++] = L[i++];
		swaps++;
		comparisons++;
	}

	/*Copy the remaining elements of R[], if there are any*/
	while (j < n2){
		array[k++] = R[j++];
		swaps++;
		comparisons++;
	}

	delete[] L, R;
}

/*******************************************************
|The heapSort function receives the array to be sorted
|and the size, the sort work by first making the array
|into a heap data structure, then it take the root of
|the tree and swap it with its last element, where now
|the last element is a sublist(the sorted array) it will
|repeat the steps and create a heap and swap it with the
|last element.
*********************************************************/
void heapSort(int array[], int size) {
	// Build heap (rearrange array) 
	for (int i = (size / 2 - 1); i >= 0; i--)
		heapify(array, size, i);

	// One by one extract an element from heap
	for (int i = size - 1; i >= 0; i--) {
		// Move current root to end
		swap(array[0], array[i]);
		swaps++;

		// call max heapify on the reduced heap
		heapify(array, i, 0);
	}
}

/*****************************************************************
|The heapify function receives an array,its size and an index
|to be converted into a heap data structure from larger to smaller
*******************************************************************/
void heapify(int array[], int size, int i) {
	int largest = i;  // Initialize largest as root
	int l = 2 * i + 1;  // left = 2*i + 1
	int r = 2 * i + 2;  // right = 2*i + 2

	// If left child is larger than root
	if (l < size && array[l] > array[largest])
	{
		largest = l;
		comparisons++;
	}

	// If right child is larger than largest so far
	if (r < size && array[r] > array[largest])
	{
		largest = r;
		comparisons++;
	}

	// If largest is not root
	if (largest != i) {
		swap(array[i], array[largest]);
		comparisons++;
		swaps++;

		// Recursively heapify the affected sub-tree
		heapify(array, size, largest);
	}
}

void quickSort(int array[], int low, int high) {
	if (low < high)
	{
		comparisons++;
		/* pi is partitioning index, arr[p] is now
		at right place */
		int pi = partition(array, low, high);

		quickSort(array, low, pi - 1);  // Before pi
		quickSort(array, pi + 1, high);  // After pi
	}
}

int partition(int array[], int low, int high) {
	// pivot (Element to be placed at right position)
	int pivot = array[high];
	int i = (low - 1);  // Index of smaller element

	for (int j = low; j <= high - 1; j++) {
		// If current element is smaller than or
		// equal to pivot
		if (array[j] <= pivot) {
			i++;    // increment index of smaller element
			swap(array[i], array[j]);
			comparisons++;
			swaps++;
		}
	}
	swap(array[i + 1], array[high]);
	swaps++;
	return (i + 1);
}

/************************************************************
|The countingSort function receives an array to be sorted,
|its size, and the range(The max element the array can have).
|This sorting algorithm work by first creating an array count
|and putting in the number of times each element appear, eg:
|if A[0] = 5 and A[3] = 5 then count will contain in position
|five,2. Then it add the previous element to the next element
|count[0] = 1, count[1] = 10 then count[1] = 11. At the end
|the sorting is a compounded function
|Sorted[count[array[i]]] = array[i] where array[i] is the index
|of count and count[array[i]] is the index of the sorted array.
***************************************************************/
void countingSort(int array[], int size, int range)
{
	//Create the index for the sorted array.
	int* count = new int[range];

	//Initialize array to 0 because there is garbage
	int i = 0;
	for (i = 0; i < range; i++)
		count[i] = 0;

	//count the number of repetition. 
	for (i = 0; i < size; i++)
		count[array[i]]++;

	//Add the previous element to the next element
	for (i = 1; i < range; i++)
		count[i] += count[i - 1];

	//create a new array to put the sorted elements in.
	int* sorted = new int[size];

	#define index count[array[i]]
	for (i = 0; i < size; i++)
	{
		sorted[index - 1] = array[i];
		index--;
	}

	//Copy sorted into array
	for (i = 0; i < size; i++)
	{
		array[i] = sorted[i];
		swaps++;
	}

	delete[] count, sorted;
}