#include <unordered_set>

using namespace std;


template<class T>
unordered_set<T> set_intersect(unordered_set<T>& a, unordered_set<T>& b)
{
	unordered_set<T> c;
	for(auto& t:a)
		if(b.count(t))
			c.insert(t);

	return c;

}