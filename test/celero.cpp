///
/// \author	John Farrier
///
/// \copyright Copyright 2014 John Farrier 
/// \copyright Copyright 2015 Dominik Köppl
///
/// Licensed under the Apache License, Version 2.0 (the "License");
/// you may not use this file except in compliance with the License.
/// You may obtain a copy of the License at
/// 
/// http://www.apache.org/licenses/LICENSE-2.0
/// 
/// Unless required by applicable law or agreed to in writing, software
/// distributed under the License is distributed on an "AS IS" BASIS,
/// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
/// See the License for the specific language governing permissions and
/// limitations under the License.
///

#include <celero/Archive.h>
#include <celero/Celero.h>
#include <celero/Console.h>
#include <celero/Benchmark.h>
#include <celero/TestVector.h>
#include <celero/Utilities.h>
#include <celero/Executor.h>
#include <celero/Print.h>
#include <celero/ResultTable.h>
#include <celero/JUnit.h>
//#include <celero/CommandLine.h>
#include <celero/Distribution.h>
#include <celero/Callbacks.h>
#include <gflags/gflags.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include "celero.hpp"

DEFINE_string(cel_group, "", "Runs a specific group of benchmarks.");
DEFINE_string(cel_outputTable, "", "Saves a results table to the named file.");
DEFINE_string(cel_junit, "", "Saves a JUnit XML-formatted file to the named file.");
DEFINE_string(cel_archive, "", "Saves or updates a result archive file.");
DEFINE_uint64(cel_distribution, 0, "Builds a file to help characterize the distribution of measurements and exits.");


void run_celero()
{
	// Initial output
	celero::print::GreenBar("");
	std::cout << std::fixed;
	celero::console::SetConsoleColor(celero::console::ConsoleColor_Green_Bold);
	std::cout << "[  CELERO  ]" << std::endl;
	celero::console::SetConsoleColor(celero::console::ConsoleColor_Default);

	celero::print::GreenBar("");
	celero::timer::CachePerformanceFrequency();
	
	// Shall we build a distribution?
	if(FLAGS_cel_distribution > 0)
	{
		celero::RunDistribution(FLAGS_cel_distribution);
	}

	// Has a result output file been specified?
	if(!FLAGS_cel_outputTable.empty())
	{
		celero::ResultTable::Instance().setFileName(FLAGS_cel_outputTable);

		celero::AddExperimentResultCompleteFunction(
			[](std::shared_ptr<celero::Result> p)
			{
				celero::ResultTable::Instance().add(p);
			});
	}

	// Has a result output file been specified?
	if(!FLAGS_cel_archive.empty())
	{
		celero::Archive::Instance().setFileName(FLAGS_cel_archive);
		
		celero::AddExperimentResultCompleteFunction(
			[](std::shared_ptr<celero::Result> p)
			{
				celero::Archive::Instance().add(p);
			});
	}

	// Has a JUnit output file been specified?
	if(!FLAGS_cel_junit.empty())
	{
		celero::JUnit::Instance().setFileName(FLAGS_cel_junit);

		celero::AddExperimentResultCompleteFunction(
			[](std::shared_ptr<celero::Result> p)
			{
				celero::JUnit::Instance().add(p);
			});
	}

	std::string finalOutput;

	// Has a run group been specified?
	if(!FLAGS_cel_group.empty())
	{
		celero::executor::Run(FLAGS_cel_group);
	}
	else
	{
		celero::executor::RunAll();
	}
	
	// Final output.
	celero::print::StageBanner(finalOutput);
}
