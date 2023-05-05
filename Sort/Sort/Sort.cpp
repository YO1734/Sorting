
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
int getRandomNumber(int min, long max)
{
	static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
	return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}

//BubbleSort
void BubbleSort(std::vector<int> &array, int size)
{
	for (int i = 0; i < size-1; i++)
	{
		for (int j = 0; j < size-1; j++)
		{
			if (array[j] > array[j+1])
			{
				std::swap(array[j], array[j +1]);
			}
			
		}

	}
}

//InsertionSort
void InsertionSort(std::vector<int>& array, int size)
{
	int key, j;
	for (int i = 1; i < size; i++)
	{
		key = array[i];
		j = i - 1;
		while (j >= 0 && array[j] > key)
		{
			std::swap(array[j + 1], array[j]);
			j--;
		}
	}
	
}

//SelectionSort
void SelectionSort(std::vector<int>& array, int size)
{
	int minIndex;
	for (int i = 0; i < size - 1; i++)
	{
		minIndex = i;
		for (int j = i + 1; j < size; j++)
		{
			if (array[j] < array[minIndex])
			{
				minIndex = j;
				
			}
		}
		if (minIndex != i)
		{
			std::swap(array[i], array[minIndex]);
		}
	}
}


//QuickSort
void QuickSort(std::vector<int>& array, int left, int right) //left- 0 index, right last index of size
{
	int i = left, j = right-1;
	int temp;
	int pivot = array[(left + right) / 2];

	/* partition */
	while (i <= j) {
		while (array[i] < pivot)
			i++;
		while (array[j] > pivot)
			j--;
		if (i <= j) {
			temp = array[i];
			array[i] = array[j];
			array[j] = temp;
			i++;
			j--;
		}
	};
	
	/* recursion */
	if (left < j)
		QuickSort(array, left, j);
	if (i < right)
		QuickSort(array, i, right);

}

//Merge Sort
void Merge(std::vector<int>& array, int left, int middle, int right) 
{
	int n1 = middle - left + 1;
	int n2 = right - middle;

	std::vector<int> leftArr(n1);
	std::vector<int> rightArr(n2);

	for (int i = 0; i < n1; i++)
		leftArr[i] = array[left + i];
	for (int j = 0; j < n2; j++)
		rightArr[j] = array[middle + 1 + j];

	int i = 0;
	int j = 0;
	int k = left;

	while (i < n1 && j < n2) {
		if (leftArr[i] <= rightArr[j]) {
			array[k] = leftArr[i];
			i++;
		}
		else {
			array[k] = rightArr[j];
			j++;
		}
		k++;
	}

	while (i < n1) {
		array[k] = leftArr[i];
		i++;
		k++;
	}

	while (j < n2) {
		array[k] = rightArr[j];
		j++;
		k++;
	}
}

void MergeSort(std::vector<int>& array, int left, int right) 
{
	if (left >= right) {
		return;
	}
	int middle = (left + right) / 2;
	MergeSort(array, left, middle);
	MergeSort(array, middle + 1, right);
	Merge(array, left, middle, right);
}

//HeapSort

void Heapify(std::vector<int>& array, int n, int i) 
{
	int largest = i;
	int l = 2 * i + 1;
	int r = 2 * i + 2;

	if (l < n && array[l] > array[largest]) {
		largest = l;
	}

	if (r < n && array[r] > array[largest]) {
		largest = r;
	}

	if (largest != i) {
		std::swap(array[i], array[largest]);
		Heapify(array, n, largest);
	}
}

void HeapSort(std::vector<int>& array) 
{
	int n = array.size();

	for (int i = n / 2 - 1; i >= 0; i--) {
		Heapify(array, n, i);
	}

	for (int i = n - 1; i >= 0; i--) {
		std::swap(array[0], array[i]);
		Heapify(array, i, 0);
	}
}

//Radix Sort
void CountSortByDigit(std::vector<int>& array, int exp) {
	std::vector<int> count(10, 0);

	// Count the occurrences of each digit
	for (int i = 0; i < array.size(); i++) {
		count[(array[i] / exp) % 10]++;
	}

	// Calculate the cumulative sum of counts
	for (int i = 1; i < 10; i++) {
		count[i] += count[i - 1];
	}

	// Build the sorted array
	std::vector<int> output(array.size());
	for (int i = array.size() - 1; i >= 0; i--) {
		output[count[(array[i] / exp) % 10] - 1] = array[i];
		count[(array[i] / exp) % 10]--;
	}

	// Copy the sorted array back to the original array
	for (int i = 0; i < array.size(); i++) {
		array[i] = output[i];
	}
}

void RadixSort(std::vector<int>& array) {
	// Determine the maximum number of digits in the largest number
	int max_digits = 0;
	for (int i = 0; i < array.size(); i++) {
		int digits = 1;
		while (array[i] / 10 > 0) {
			digits++;
			array[i] /= 10;
		}
		max_digits = std::max(max_digits, digits);
	}

	// Sort the array by each digit, starting with the least significant
	for (int exp = 1; exp <= std::pow(10, max_digits - 1); exp *= 10) {
		CountSortByDigit(array, exp);
	}
}

//Exchange Sort
void ExchangeSort(std::vector<int>& array)
{
	for (int i = 0; i < array.size(); i++)
	{
		for (int j = i + 1; j < array.size(); j++)
		{
			if (array[i] > array[j])
			{
				std::swap(array[i], array[j]);
			}
		}
	}
}

//ShellSort

void ShellSort(std::vector<int>& array)
{
	int n = array.size();
	for (int gap = n / 2; gap > 0; gap /= 2)
	{
		for (int i = gap; i < n; i++)
		{
			int temp = array[i];
			int j;
			for (j = i; j >= gap && array[j - gap] > temp; j -= gap)
			{
				array[j] = array[j - gap];
			}
			array[j] = temp;
		}
	}
}

//BucketSort

void BucketSort(std::vector<int>& array) {
	int n = array.size();

	// Создаем массив корзин
	std::vector<std::vector<int>> buckets(n);

	// Размещаем элементы в соответствующие корзины
	for (int i = 0; i < n; ++i) {
		int bucket_idx = n * (static_cast<double>(array[i]) / (std::numeric_limits<int>::max() + 1.0));
		buckets[bucket_idx].push_back(array[i]);
	}

	// Сортируем каждую корзину
	for (int i = 0; i < n; ++i) {
		sort(buckets[i].begin(), buckets[i].end());
	}

	// Склеиваем корзины вместе
	int idx = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < buckets[i].size(); ++j) {
			array[idx++] = buckets[i][j];
		}
	}
}

void FillArray(std::vector<int> &array, int size)
{
	for (int i = 0; i < size; i++)
	{
		
		array[i] = getRandomNumber(0,100) ;
	}
}


void PrintArray(std::vector<int> &array, int size)
{
	for (int i=0;i<size;i++)
	{
		std::cout << array[i] << std::endl;
	}

}

int main()
{

	srand(static_cast<unsigned int>(time(0)));
	rand();


	std::vector <int>myarray;
	myarray.resize(20);

	FillArray(myarray, myarray.size());
	//BubbleSort(myarray,myarray.size());
	//InsertionSort(myarray, myarray.size());
	//QuickSort(myarray, 0, myarray.size());
	//SelectionSort(myarray, myarray.size());
	//MergeSort(myarray, 0, myarray.size() - 1);
	//HeapSort(myarray);
	//RadixSort(myarray);
	//ExchangeSort(myarray);
	//ShellSort(myarray);
	BucketSort(myarray);
	PrintArray(myarray, myarray.size());


}

