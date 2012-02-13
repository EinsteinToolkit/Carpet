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

#include <assert.h>
#include <map>
#include <string>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include "TimerNode.hh"

namespace Carpet {

  using namespace std;

  TimerNode *TimerNode::d_current = 0;
  TimerNode *TimerNode::root_timer = 0;

  TimerNode::TimerNode(string name): d_name(name), d_parent(0),
                                     d_running(false)
  {
    timer = 0;
  }

  TimerNode::~TimerNode()
  {
    for(map<string,TimerNode*>::iterator iter=d_children.begin();
        iter!=d_children.end(); ++iter)
    {
      delete iter->second;
    }
    d_children.clear();
    delete timer;
  }

  string TimerNode::pathName() const
  {
    if (d_parent)
      return d_parent->pathName() + string("/") + getName();
    else
      return getName();
  }

  void TimerNode::start()
  {
    if(!d_running)
    {
      d_running = true;
      d_parent = d_current;
      d_current=this;
      if (timer == 0)
        timer = new CactusTimer(pathName());
      assert(timer);
      timer->start();
    }
    else
    {
      assert(0); // Timer is already running
    }
  }

  void TimerNode::stop()
  {
    if(d_running)
    {
      // A timer can only be stopped if it is the current timer
      if(this != d_current)
        CCTK_VWarn (0, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Tried to stop non-current timer '%s'", getName().c_str());

      timer->stop();

      d_running=false;
      d_current=d_parent;
    }
    else
    {
      assert(0); // Timer is not running
    }
  }

  /// Get the name of the timer
  string TimerNode::getName() const
  {
    assert(d_name.length() > 0);
    return d_name;
  }

  /// Get the current timer
  TimerNode* TimerNode::getCurrentTimer()
  {
    return d_current;
  }

  /// Get the root of the timer tree
  TimerNode* TimerNode::getRootTimer()
  {
    if (root_timer == 0)
      root_timer = new TimerNode("main");
    return root_timer;
  }

  /// Find the child timer that matches the name provided.  If it is
  /// not found then a new timer with that name is allocated.
  TimerNode* TimerNode::getChildTimer(string name)
  {
    // Find child
    TimerNode *child=d_children[name];

    // If the pointer is null then allocate it
    if(child == 0)
      d_children[name] = child = new TimerNode(name);

    return child;
  }

  /// Get the time measured by this timer
  double TimerNode::getTime()
  {
    return timer->getTime();
  }

  /// Print this node and its children as an ASCII tree
  void TimerNode::print(ostream& out, double total, int level,
                        double threshold)
  {
    string space;

    // Compute the level of indentation for this depth
    for(int i=0;i<level-1;i++)
      space += "| ";

    if (level != 0)
      space += "|_";

    const double t = getTime();

    if (level == 0)
    {
      out << "-----------------------" << endl;
      out << "Percent   t/secs  Timer" << endl;
      out << "-----------------------" << endl;
    }

    const int pcw = 6;
    const int tw = 8;

    // Print this timer value
    out << fixed << setw(pcw) << setprecision(1) << 100.0 * t / total << "%"
        << " " << fixed << setw(tw) << setprecision(1) << t
        << "  " << space << d_name << endl;

    double children_time = 0;
    bool printed_children = false;

    // Recursively print the children
    for(map<string,TimerNode*>::iterator iter = d_children.begin();
        iter != d_children.end(); iter++)
    {
      if (iter->second->getTime() * 100.0 / total > threshold)
      {
        iter->second->print(out,total,level+1,threshold);
        printed_children = true;
      }
      children_time += iter->second->getTime();
    }

    if (d_children.size() > 0 && printed_children)
    {
      const double untimed = t - children_time;
      
      if (100.0 * untimed / total > threshold)
      {
        // Print the untimed portion
        out << fixed << setw(pcw) << setprecision(1) << 100.0 * untimed / total << "%"
            << " " << fixed << setw(tw) << setprecision(1) << untimed
            << "  | " << space << "untimed" << endl;
      }
    }
  }

  /// Print this node and its children as an XML file
  void TimerNode::printXML(ostream& out, int level)
  {
    string space;

    // Compute the level of indentation for this node
    for(int i=0;i<level;i++)
      space=space+"  ";

    out << space << "<timer name = " << "\"" << escapeForXML(d_name) << "\"> ";
    out << getTime() << " ";

    // For compactness, only use multiple lines if there are children
    if (d_children.size() != 0)
    {
      out << endl;

      // Recursively print the children
      for(map<string,TimerNode*>::iterator iter=d_children.begin();iter!=d_children.end();++iter)
        iter->second->printXML(out,level+1);
      out << space;
    }

    out << "</timer>" << endl;
  }

  /// Make a string suitable for inclusion in an XML file
  string TimerNode::escapeForXML(const string &s) const
  {
    // XML attributes cannot contain unescaped angle-brackets.  As a
    // simple solution to this, replace them with | characters.

    string s2(s);
    using std::string;
    using std::cout;
    using std::endl;
    using std::replace;

    replace(s2.begin(), s2.end(), '<', '|');
    replace(s2.begin(), s2.end(), '>', '|');
    replace(s2.begin(), s2.end(), '&', '|');

    return s2;
  }
}
