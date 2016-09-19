//
// Created by Bernardo Clavijo (TGAC) on 19/09/2016.
//

#ifndef W2RAP_CONTIGGER_W2RAP_TIMERS_H
#define W2RAP_CONTIGGER_W2RAP_TIMERS_H

#ifdef TIME_LOGGING
#include <string>
#include <sstream>
#include <vector>
#include <chrono>


//sets the timer variables
#define TIMELOG_DECLARE(NAME) uint64_t NAME=0,NAME ## starttime=0;

#define TIMELOG_DECLARE_LOCAL(NAME,SUFFIX) uint64_t NAME ## SUFFIX ## starttime=0;

#define TIMELOG_DECLARE_ATOMIC(NAME) std::atomic_uint_fast64_t NAME(0),NAME ## starttime(0);

#define TIMELOG_CREATE_GLOBAL(NAME) std::atomic_uint_fast64_t NAME(0);
#define TIMELOG_DECLARE_GLOBAL(NAME) extern std::atomic_uint_fast64_t NAME;

//saves the start of a timer (in TIMER_(TIMER_NAME)_start variable)
#define TIMELOG_START(NAME) {NAME ## starttime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();}

#define TIMELOG_START_LOCAL(NAME,SUFFIX) {NAME ## SUFFIX ## starttime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();}


//computes the time since the start, and adds it to the timer
#define TIMELOG_STOP(NAME) {NAME+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count()-NAME ## starttime;}

#define TIMELOG_STOP_LOCAL(NAME,SUFFIX) {NAME+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count()-NAME ## SUFFIX ## starttime;}

//adds a partial time to a global timer
#define TIMELOG_ADD_TO(DEST, SOURCE) {DEST+=SOURCE;}

#define TIMELOG_RESET(NAME) {NAME=0;}

//shows the total of each timer, and its % from the total of timers
#define TIMELOG_REPORT(OSDEST,TITLE,...) { std::stringstream ss(#__VA_ARGS__);\
std::string item;\
uint64_t timers[] = { __VA_ARGS__ };\
uint64_t total = 0; for (auto &t:timers) total+=t;\
uint64_t timers_perc[] = { __VA_ARGS__ };\
for (auto &t: timers_perc) t=t*100/total;\
OSDEST<<"TIME REPORT FOR '"<<#TITLE<<"': Total time: "<<total<<std::endl;\
int i=0;\
while (std::getline(ss, item, ',')) {OSDEST<<item<<": "<<timers[i]<<" ("<<timers_perc[i]<<"%)   "; ++i;}\
OSDEST<<std::endl;}

#endif

#ifndef TIME_LOGGING
#define TIMELOG_DECLARE(NAME)
#define TIMELOG_DECLARE_LOCAL(NAME,SUFFIX)
#define TIMELOG_DECLARE_ATOMIC(NAME)
#define TIMELOG_CREATE_GLOBAL(NAME)
#define TIMELOG_DECLARE_GLOBAL(NAME)
#define TIMELOG_START(NAME)
#define TIMELOG_START_LOCAL(NAME,SUFFIX)
#define TIMELOG_STOP(NAME)
#define TIMELOG_STOP_LOCAL(NAME,SUFFIX)
#define TIMELOG_ADD_TO(DEST, SOURCE) {DEST+=SOURCE;}
#define TIMELOG_RESET(NAME)
#define TIMELOG_REPORT(OSDEST,TITLE,...)
#endif

//********** GLOBAL TIMERS DECLARATION ************
// Only declare timers here if they need to be global, ideally only if you intend to use them across threads.
TIMELOG_DECLARE_GLOBAL(CP1_Align);
TIMELOG_DECLARE_GLOBAL(CP1_MakeStacks);
TIMELOG_DECLARE_GLOBAL(CP1_Correct);
TIMELOG_DECLARE_GLOBAL(C1P_Align);
TIMELOG_DECLARE_GLOBAL(C1P_InitBasesQuals);
TIMELOG_DECLARE_GLOBAL(C1P_Correct);
TIMELOG_DECLARE_GLOBAL(C1P_UpdateBasesQuals);
#endif //W2RAP_CONTIGGER_W2RAP_TIMERS_H
