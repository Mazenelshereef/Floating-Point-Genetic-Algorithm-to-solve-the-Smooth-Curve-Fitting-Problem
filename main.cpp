#include <cmath>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <bits/stdc++.h>
#include <utility>
#include <iostream>
#include <set>

using namespace std;
typedef long long unsigned int ll;
ll numtest;
float calfitness(vector<float> chromosome, vector<pair<float, float>> vec)
{
    float sum = 0;
    float counter = 0;

    for (size_t i = 0; i < vec.size(); i++)
    {
        counter = 0;
        for (size_t j = 0; j < chromosome.size(); j++)
        {
            counter += chromosome[j] * powf(vec[i].first, j);
        }
        counter -= vec[i].second;
        sum += powf(counter, 2);
    }
    return 1.0 / (sum / vec.size());
}
vector<float> mutation(vector<float> chromosome, int generation)
{
    for (size_t i = 0; i < chromosome.size(); i++)
    {
        float lx, ux, r1, y, delta, r;
        float pm = 0.1;
        srand(static_cast<unsigned int>(time(nullptr)));
        float randNum = rand() % (1 - 0 + 1) + 0;
        if (randNum <= pm)
        {
            lx = chromosome[i] + 10;
            ux = 10 - chromosome[i];
            r1 = (float)rand() / RAND_MAX;
            if (r1 > 0.5)
                y = lx;
            else
                y = ux;
            r = rand() % (1 - 0 + 1) + 0;
            delta = y * (1 - powf(r, (1 - generation / (0.01))));
            if (r1 < 0.5)
                chromosome[i] = chromosome[i] - delta;
            else
                chromosome[i] = chromosome[i] + delta;
        }
    }
    return chromosome;
}
pair<float, vector<float>> tournment(vector<pair<float, vector<float>>> pop, int k)
{
    pair<float, vector<float>> bestfit = pop[1];
    set<int> cand;
    srand(static_cast<unsigned int>(time(nullptr)));
    while (cand.size() != k)
    {
        int randNum = rand() % (pop.size() - 0 + 1) + 0;
        cand.insert(randNum);
    }
    for (size_t i = 1; i < cand.size(); i++)
    {
        if (pop[i].first > bestfit.first)
        {
            bestfit = pop[i];
        }
    }
    return bestfit;
}

ll numOfData, pynomialDegree, numItems;
const int populationSize = 100;
float x, y, counter = 0;
float random, parent1, parent2, numGenerations = 5;
pair<float, vector<float>> firstchild1, secondchild2;
vector<pair<float, vector<float>>> candidates;
vector<float> chromosomesFitness, chromosomesCumulativeFitness;
vector<float> firstparent, secondparent, firstchild, secondchild;
vector<std::pair<float, float>> vec;
int main()
{
    freopen("input.txt","r",stdin);
    freopen("output.txt","w",stdout);

    cin >> numtest;
    for (int z = 0; z < numtest; ++z)
    {
        counter = 0;
        vector<pair<float, vector<float>>> chromosome(populationSize);
        firstparent.clear();
        secondparent.clear();
        firstchild.clear();
        secondchild.clear();
        candidates.clear();
        chromosomesFitness.clear();
        chromosomesCumulativeFitness.clear();
        vec.clear();
        pair<float, float> pair;
        cin >> numOfData;
        cin >> pynomialDegree;
        numItems = pynomialDegree + 1;
        for (int i = 0; i < numOfData; ++i)
        {
            cin >> x >> y;
            pair = make_pair(x, y);
            vec.push_back(pair);
        }
        srand(static_cast<unsigned int>(time(nullptr)));
        for (int k = 0; k < populationSize; k++)
        {

            for (int j = 0; j < numItems; j++)
            {
                float random = -10 + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (10 - (-10)))); // random  0 or 1

                chromosome[k].second.push_back(random);
            }
        }
        // applying fitness function
        for (int j = 0; j < populationSize; j++)
        {
            counter = 0;
            counter = calfitness(chromosome[j].second, vec);
            chromosome[j].first = counter;
            chromosomesFitness.push_back(counter);
        }
        // loop for generations
        for (int g = 0; g < numGenerations; g++)
        {
            firstparent.clear(); secondparent.clear();
            firstparent = tournment(chromosome, 8).second;
            secondparent = tournment(chromosome, 8).second;

            srand(static_cast<unsigned int>(time(nullptr)));
            int xc = rand() % (numItems - 0 + 1) + 0;
            int dc = xc;
            while (dc == xc)
            { // make sure dc not equal xc
                dc = rand() % (numItems - 0 + 1) + 0;
            }
            if (xc > dc)
            { // make sure dc always bigger than xc
                int temp = xc;
                xc = dc;
                dc = temp;
            }
            int rc = rand() % (numItems - 2 + 1) + 1; // random between 0 and 1
            int pc = 0.7;
            int pm = 0.1;

            if (rc <= pc) // crossover
            {
                for (int i = 0; i < xc; i++)
                {
                    firstchild[i] = firstparent[i];
                    secondchild[i] = secondparent[i];
                }
                for (int i = xc; i < dc; i++)
                {
                    firstchild[i] = secondparent[i];
                    secondchild[i] = firstparent[i];
                }
                for (int i = dc; i < firstparent.size(); i++)
                {
                    firstchild[i] = firstparent[i];
                    secondchild[i] = secondparent[i];
                }
            }
            else // no crossover
            {
                firstchild = firstparent;
                secondchild = secondparent;
            }

            firstchild = mutation(firstchild, g);
            firstchild = mutation(secondchild, g);
            std::sort(chromosome.begin(), chromosome.end());
            float ff = calfitness(firstchild, vec);  // first fitness
            float sf = calfitness(secondchild, vec); // second fitness
            firstchild1.first = ff;
            secondchild2.first = sf;
            firstchild1.second = firstchild;
            secondchild2.second = secondchild;

            candidates.push_back(firstchild1);
            candidates.push_back(secondchild2);
            candidates.push_back(chromosome[0]);
            candidates.push_back(chromosome[1]);
            sort(candidates.begin(), candidates.end());
            if (true)
            {
                chromosome[0] = candidates[2];
            }
            if (true)
            {
                chromosome[1] = candidates[3];
            }
        }
        std::sort(chromosome.begin(), chromosome.end());
        vector<float> bestOne = chromosome[populationSize - 1].second;
        float finalFitness = calfitness(bestOne, vec);

        cout << "" << endl;
        cout << "---------------------"
             << "TestCase " << z + 1 << "---------------------" << endl;
        for (size_t i = 0; i < bestOne.size(); i++)
        {
            cout << "a[" << i << "]"
                 << " = " << bestOne[i] << endl;
        }

        cout << "MSE = " << 1 / calfitness(bestOne, vec);
    }
    return 0;
}
