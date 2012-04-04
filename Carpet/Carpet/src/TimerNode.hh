/*
 *
 * The MIT License
 *
 * Copyright (c) 1997-2010 Center for the Simulation of Accidental Fires and
 * Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI),
 * University of Utah.
 *
 * License for the specific language governing rights and limitations under
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * */

/* This class was originally written by Justin Luitjens and
   subsequently integrated with Cactus/Carpet by Ian Hinder and
   heavily modified. */

#ifndef TIMERNODE_HH
#define TIMERNODE_HH

#include <assert.h>
#include <map>
#include <string>
#include <ostream>
#include <iostream>

#include "CactusTimer.hh"

namespace Carpet {

using namespace std;

/**

The TimerNode class implements a tree structure where each node
represents a timer, implemented as a CactusTimer.  Each node of the
tree can have zero of more children, where the names of the child
nodes are unique within a single parent, but not necessarily unique
within the entire tree.  A tree formed of TimerNode objects represents
the execution profile of the program, inasmuch as it is instrumented
by timers.

Child nodes of a given name are accessed using the getChildTimer
method, where a node is created with the given name if none exists
already.  This ensures that the names of the child timers are unique.

*/

class TimerNode
{
public:
  TimerNode(string name);
  ~TimerNode();

  void start();
  void stop();

  string getName() const;
  string pathName() const;

  static TimerNode* getCurrentTimer();

  // Finds the child timer that matches the name provided.  If it is
  // not found then that timer is allocated.
  TimerNode* getChildTimer(string name);
  static TimerNode* getRootTimer();

  double getTime();

  void print(ostream& out, double total, int level=0, double threshold=0.0, int precision=1);
  void printXML(ostream& out, int level=0);

private:
  string escapeForXML(const string &s) const;

  string d_name;
  std::map<string,TimerNode*> d_children;
  static TimerNode* d_current;
  static TimerNode* root_timer;
  TimerNode *d_parent;
  bool d_running;
  CactusTimer *timer;
};
}

#endif // TIMERNODE_HH
