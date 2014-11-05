#include <unordered_set>
#include <algorithm>
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


template <class T> 
vector<T> set2Vector(unordered_set<T>& a)
{
	int count =0;

	vector<T> result(a.size());
	for(auto&s:a)	
		result[count++]=s;

	sort(result.begin(),result.end());

	return result;
}


template<class T>
unordered_set<T> vector2Set(const vector<T>& a)
{
	unordered_set<T> result;
	for (auto& s:a)
	{
		result.insert(s);
	}

	return result;
}
