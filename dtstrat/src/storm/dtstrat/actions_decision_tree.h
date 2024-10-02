
/*
This inline function was automatically generated using DecisionTreeToCpp Converter
It takes feature vector as single argument:
feature_vector[0] - x0
feature_vector[1] - x1

It returns the predicted class:
0 - 0.0
1 - 1.0
2 - 3.0
3 - 5.0
4 - 7.0
5 - 10.0
6 - 11.0
7 - 12.0
8 - 13.0
9 - 14.0
10 - 16.0

Simply include this file to your project and use it
*/
#include <vector>

inline int actions_decision_tree(const std::vector<double> & feature_vector) 
{
	if (feature_vector.at(1) <= 1.5) {
		if (feature_vector.at(1) <= 0.5) {
			if (feature_vector.at(0) <= 0.5) {
				return 11.0;
			}
			else {
				if (feature_vector.at(0) <= 1.5) {
					return 10.0;
				}
				else {
					return 11.0;
				}
			}
		}
		else {
			if (feature_vector.at(0) <= 0.5) {
				return 1.0;
			}
			else {
				return 12.0;
			}
		}
	}
	else {
		if (feature_vector.at(0) <= 1.5) {
			if (feature_vector.at(0) <= 0.5) {
				return 13.0;
			}
			else {
				return 3.0;
			}
		}
		else {
			if (feature_vector.at(1) <= 4.5) {
				if (feature_vector.at(1) <= 3.5) {
					if (feature_vector.at(1) <= 2.5) {
						if (feature_vector.at(0) <= 2.5) {
							return 0.0;
						}
						else {
							if (feature_vector.at(0) <= 3.5) {
								return 5.0;
							}
							else {
								if (feature_vector.at(0) <= 4.5) {
									return 7.0;
								}
								else {
									return 0.0;
								}
							}
						}
					}
					else {
						return 14.0;
					}
				}
				else {
					return 16.0;
				}
			}
			else {
				return 0.0;
			}
		}
	}
}