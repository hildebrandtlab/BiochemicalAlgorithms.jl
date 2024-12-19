// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: AmberFF_bench.C,v 1.5 2002/06/16 19:12:55 oliver Exp $
#include <BALLBenchmarkConfig.h>
#include <BALL/CONCEPT/benchmark.h>

///////////////////////////

#include <BALL/MOLMEC/AMBER/amber.h>
#include <BALL/MOLMEC/COMMON/forceFieldComponent.h>
#include <BALL/FORMAT/PDBFile.h>

///////////////////////////

using namespace BALL;

START_BENCHMARK(AmberFF, 1.0, "$Id: AmberFF_bench.C,v 1.5 2002/06/16 19:12:55 oliver Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PDBFile pdb(BALL_BENCHMARK_DATA_PATH(AmberFF_bench.pdb));
System S;
pdb >> S;

int n = 1000;

AmberFF amber;
START_SECTION(1000x setup, 1.0)
	START_TIMER
		for (int i = 0; i < n; i++)
		{
			amber.setup(S);
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x update w/o selection, 1.0)
	START_TIMER
		for (int i = 0; i < n; i++)
		{
			amber.update();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x energy calculation w/o selection, 1.0)
	START_TIMER
		for (int i = 0; i < n; i++)
		{
			amber.updateEnergy();	
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x force calculation w/o selection, 1.0)
	START_TIMER
		for (int i = 0; i < n; i++)
		{
			amber.updateForces();
		}
	STOP_TIMER
END_SECTION

ForceFieldComponent* component;
START_SECTION(1000x nonbonded energy calculation w/o selection, 1.0)
	component = amber.getComponent("Amber NonBonded");
	START_TIMER
		for (Size i = 0; i < n; i++)
		{
			component->updateEnergy();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x stretch energy calculation w/o selection, 1.0)
	component = amber.getComponent("Amber Stretch");
	START_TIMER
		for (Size i = 0; i < n; i++)
		{
			component->updateEnergy();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x bend energy calculation w/o selection, 1.0)
	component = amber.getComponent("Amber Bend");
	START_TIMER
		for (Size i = 0; i < n; i++)
		{
			component->updateEnergy();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x torsion energy calculation w/o selection, 1.0)
	component = amber.getComponent("Amber Torsion");
	START_TIMER
		for (Size i = 0; i < n; i++)
		{
			component->updateEnergy();
		}
	STOP_TIMER
END_SECTION

S.beginResidue()->select();
START_SECTION(1000x update w/ selection, 1.0)
	START_TIMER
		for (int i = 0; i < n; i++)
		{
			amber.update();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x energy calculation w/ selection, 1.0)
	START_TIMER
		for (int i = 0; i < n; i++)
		{
			amber.updateEnergy();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x force calculation w/ selection, 1.0)
	START_TIMER
		for (int i = 0; i < n; i++)
		{
			amber.updateForces();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x nonbonded energy calculation w/ selection, 1.0)
	component = amber.getComponent("Amber NonBonded");
	START_TIMER
		for (Size i = 0; i < n; i++)
		{
			component->updateEnergy();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x stretch energy calculation w/ selection, 1.0)
	component = amber.getComponent("Amber Stretch");
	START_TIMER
		for (Size i = 0; i < n; i++)
		{
			component->updateEnergy();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x bend energy calculation w/ selection, 1.0)
	component = amber.getComponent("Amber Bend");
	START_TIMER
		for (Size i = 0; i < n; i++)
		{
			component->updateEnergy();
		}
	STOP_TIMER
END_SECTION

START_SECTION(1000x torsion energy calculation w/ selection, 1.0)
	component = amber.getComponent("Amber Torsion");
	START_TIMER
		for (Size i = 0; i < n; i++)
		{
			component->updateEnergy();
		}
	STOP_TIMER
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_BENCHMARK
