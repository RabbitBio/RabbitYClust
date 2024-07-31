#include "KHFMinHash.h"

#include <iostream>
#include <string>

using namespace std;
using namespace Sketch;

int main()
{
	char test[] = "LKPYLDKGDIIIDGGNTFFQDTIRRNRELSAEGFNFIGTGVSGGEEGALKGPSIMPGGQKEAYELVAPILTKIAAVAEDGLKPYLDKGDIIIDGGNTFFQDTIRRNRELSAEGFNFIGTGVSGGEEGALKGPSIMPGGQKEAYELVAPILTKIAAVAEDG";

	KHFMinHash mh = KHFMinHash(test);			
	
	auto & sketch = mh.getSektch();	

	for(int i = 0; i < sketch.hashes.size(); i++)
		cout << sketch.hashes[i] << endl;

	return 0;
	
}
