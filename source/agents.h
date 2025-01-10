#pragma once
#include "helpers.h"
#include "grid_agent.forward.h"
#include "kernel.h"


const float AGB_coeff_a = -0.0299f;
const float AGB_coeff_b = 2.673f;

using namespace help;

class Strategy {
public:
	Strategy() = default;
	Strategy(string _vector, float _seed_mass, float _diaspore_mass, int _no_seeds_per_diaspore, float _seed_tspeed,
		float _pulp_to_seed_ratio, float _recruitment_probability, float _seedling_dbh, float _relative_growth_rate, float _seed_reserve_mass
	) :
		vector(_vector), seed_mass(_seed_mass), diaspore_mass(_diaspore_mass), no_seeds_per_diaspore(_no_seeds_per_diaspore),
		seed_tspeed(_seed_tspeed), pulp_to_seed_ratio(_pulp_to_seed_ratio), recruitment_probability(_recruitment_probability),
		seedling_dbh(_seedling_dbh), relative_growth_rate(_relative_growth_rate), seed_reserve_mass(_seed_reserve_mass)
	{}
	void print() {
		printf("id: %d, seed_mass: %f, diaspore_mass: %f, no_seeds_per_diaspore: %d, vector: %s, pulp to seed ratio: %f, seed terminal speed: %f, germination prob: %f\n",
			id, seed_mass, diaspore_mass, no_seeds_per_diaspore, vector.c_str(), pulp_to_seed_ratio, seed_tspeed, recruitment_probability
		);
	}
	int id = -1;
	float seed_mass = 0;
	float diaspore_mass = 0;
	int no_seeds_per_diaspore = 0;
	float seed_tspeed = 0;
	float pulp_to_seed_ratio = 0;
	float recruitment_probability = 0;
	float relative_growth_rate = 0;
	float seed_reserve_mass = 0;
	float seedling_dbh = 0;
	string vector = "none";
};


class StrategyGenerator {
public:
	StrategyGenerator() = default;
	StrategyGenerator(map<string, map<string, float>> user_parameters) {
		for (auto& [trait, value_range] : user_parameters) {
			ProbModel prob_model;
			string distribution_type = distribution_types[value_range["distribution_type"]];
			if (distribution_type == "uniform") {
				prob_model = ProbModel(value_range["min"], value_range["max"]);
			}
			else if (distribution_type == "linear") {
				prob_model = ProbModel(value_range["q1"], value_range["q2"], value_range["min"], value_range["max"]);
			}
			else if (distribution_type == "normal") {
				prob_model = ProbModel(value_range["mean"], value_range["stdev"], value_range["min"], value_range["max"], 0);
			}
			else if (distribution_type == "discrete") {
				prob_model = ProbModel(value_range["probability0"], value_range["probability1"], value_range["probability2"]);
			}
			else if (distribution_type == "constant") {
				prob_model = ProbModel(value_range["constant"]);
			}
			trait_distributions[trait] = prob_model;
		}
	}
	string pick_vector() {
		int vector = trait_distributions["vector"].sample();
		if (vector == 0) return "linear";
		else if (vector == 1) return "wind";
		else return "animal";
	}
	float compute_wing_mass(float cumulative_seed_mass) {
		// Estimate wing mass based on correlation we found in seed- and wing measurement data taken from (Greene and Johnson, 1993)
		return max(0, (cumulative_seed_mass - 0.03387f) / 3.6609f);
	}
	float compute_diaspore_mass(float no_seeds_per_diaspore, float seed_mass, string vector, float fruit_pulp_mass, float &pulp_to_seed_ratio) {
		float cumulative_seed_mass = seed_mass * no_seeds_per_diaspore;
		float diaspore_mass;
		if (vector == "wind") {
			float wing_mass = compute_wing_mass(cumulative_seed_mass);
			diaspore_mass = cumulative_seed_mass + wing_mass;
			//printf("cumulative seed mass: %f, wing mass: %f \n", cumulative_seed_mass, wing_mass);
		}
		else if (vector == "animal") {
			diaspore_mass = cumulative_seed_mass + fruit_pulp_mass;
		}
		else {
			diaspore_mass = cumulative_seed_mass;
		}
		pulp_to_seed_ratio = fruit_pulp_mass / cumulative_seed_mass;
		return diaspore_mass;
	}
	float calculate_tspeed(float diaspore_mass) {
		// Compute terminal descent velocity based on correlation presented by (Greene and Johnson, 1993)
		return 0.501f * pow(diaspore_mass * 1000, 0.174);
	}
	int sample_no_seeds_per_diaspore() {
		int no_seeds_per_diaspore = round(trait_distributions["no_seeds_per_diaspore"].sample());
		if (no_seeds_per_diaspore < 1) no_seeds_per_diaspore = 1; // Ensure at least one seed per diaspore
		return no_seeds_per_diaspore;
	}
	float sample_seed_mass(string vector) {
		float seed_mass = 1.0f;
		if (vector == "wind") {
			seed_mass = trait_distributions["seed_mass_wind"].sample();
		}
		if (vector == "animal") seed_mass = trait_distributions["seed_mass_animal"].sample();
		return seed_mass;
	}
	float sample_fruit_pulp_mass() {
		return trait_distributions["fruit_pulp_mass"].sample();
	}
	float calculate_recruitment_probability(float seed_mass) {
		float proportion_eaten = 1.0f / (1.0f + exp(2.49f - 13.0f * seed_mass)); // Relation between seed mass and proportion of seeds eaten by rodents, from Wang and Ives (2017), figure 5
		float seedling_establishment_probability = max(0, 0.034 * exp(3.3 * seed_mass)); // Fitted to data from Barczyk et al (2024), see file 'seed weight vs seedling success.xlsx'
		return (1.0f - proportion_eaten) * seedling_establishment_probability;
	}
	float get_seed_reserve_mass(float seed_mass) {
		float seed_reserve_mass = 0.7313f * pow(seed_mass, 1.0633f);	// In grams, fitted to Boot (1994) data, table 2a (see file 'Relationship seed mass to seedling diameter.xlsx').
		return seed_reserve_mass;
	}
	float get_relative_growth_rate(float seed_mass) {
		float relative_growth_rate = (1.0f - 6.97f * log10(seed_mass)) * 0.001f;		// Relative growth rate of seedlings in g/g/day, based on Rose (2003), 
																						// table 2.2 (regression slope 'overall'). We multiply by 0.001 to convert mg/g/day to g/g/day.
		return relative_growth_rate;
	}
	float calculate_seedling_dbh(float seed_reserve_mass, float relative_growth_rate, float seed_mass, string vector) {
		if (trait_distributions["seedling_dbh_" + vector].sample() >= 0) {
			return trait_distributions["seedling_dbh_" + vector].sample(); // If a constant value is provided in the parameter file, use that value.
		}

		float new_mass_grams = pow(10.0f, seed_reserve_mass * exp(relative_growth_rate * 365.25f));				// Mass in grams after 1 year of growth, based on Rose (2003). 
																													// Seedling start mass is assumed to be equal to seed reserve mass.
		
		float Leaf_Mass_Fraction = 0.4f - 0.05f * log10(seed_mass);												// g/g. Based on table 2.2, Rose (2003).
		float Leaf_Area_Ratio = 10.0f - 4.93f * log10(seed_mass);												// m^2 / kg. Based on table 2.2, Rose (2003).
		float density_factor = 1.0f / (Leaf_Area_Ratio / 6.87f);												// Compute factor to correct for differences in mass density.
																													// We assume an average seed mass of 4.31 g (Figure 2.1, Rose (2003)) 
																													// and thus normalize LAR to LAR_norm = LAR_{seedmass=4.31}.
																													// LAR_norm will increase for smaller seeds (more m^2 leaf area per kg)
																													// and decrease for larger seeds. The density factor varies inversely.
																												
		//printf("seed_reserve_mass: %f, seed mass: %f, density factor: %f \n", seed_reserve_mass, seed_mass, density_factor);
		new_mass_grams /= density_factor;																		// Correct for differences in mass density (rapidly growing seedlings have
																													// lower density and should thus have a larger dbh for the same mass).												
		float new_AGB_kilograms = 0.000666f * new_mass_grams;													// Aboveground mass in kilograms. We assume a shoot-to-root ratio of 2/3
		float AGB_c = -(log(new_AGB_kilograms) + 2.5f);															// Derived from model 7, Chave et al (2014)
		float ln_new_dbh = help::get_lowest_solution_for_quadratic(0, AGB_coeff_a, AGB_coeff_b, AGB_c);			// Solve for ln(dbh) using quadratic formula (derived from model 7, 
																													// Chave et al (2014))										

		float firstyear_dbh = exp(ln_new_dbh);
		//printf("first year dbh: %f, new mass AGB: %f, relative growth rate (g/g/day): %f \n", firstyear_dbh, new_AGB_kilograms, relative_growth_rate);

		return firstyear_dbh;
	}

	void generate(Strategy &strategy) {
		int no_seeds_per_diaspore = sample_no_seeds_per_diaspore();
		string vector = pick_vector();
		float seed_mass = sample_seed_mass(vector);
		float fruit_pulp_mass = sample_fruit_pulp_mass();
		float pulp_to_seed_ratio;
		float diaspore_mass = compute_diaspore_mass(no_seeds_per_diaspore, seed_mass, vector, fruit_pulp_mass, pulp_to_seed_ratio);
		float seed_tspeed = calculate_tspeed(diaspore_mass);
		float recruitment_probability = calculate_recruitment_probability(seed_mass);
		float seed_reserve_mass = get_seed_reserve_mass(seed_mass);
		float relative_growth_rate = get_relative_growth_rate(seed_mass);
		float seedling_dbh = calculate_seedling_dbh(seed_reserve_mass, relative_growth_rate, seed_mass, vector);
		strategy = Strategy(
			vector, seed_mass, diaspore_mass, no_seeds_per_diaspore, seed_tspeed, pulp_to_seed_ratio, recruitment_probability,
			seedling_dbh, relative_growth_rate, seed_reserve_mass
		);
	}
	void mutate(Strategy& strategy, float mutation_rate) {
		bool do_mutation = help::get_rand_float(0, 1) < mutation_rate;
		if (do_mutation) {
			int trait_idx = help::get_rand_int(0, 3);
			float fruit_pulp_mass = 0;
			if (trait_idx == 0) {
				strategy.seed_mass = sample_seed_mass(strategy.vector);
			}
			else if (trait_idx == 1) {
				strategy.no_seeds_per_diaspore = sample_no_seeds_per_diaspore();
			}
			else if (trait_idx == 2) {
				fruit_pulp_mass = sample_fruit_pulp_mass();
			}
			else if (trait_idx == 3) {
				strategy.vector = pick_vector();
			}
			float pulp_to_seed_ratio;
			strategy.diaspore_mass = compute_diaspore_mass(strategy.no_seeds_per_diaspore, strategy.seed_mass, strategy.vector, fruit_pulp_mass, pulp_to_seed_ratio);
			strategy.seed_tspeed = calculate_tspeed(strategy.diaspore_mass);
			strategy.pulp_to_seed_ratio = pulp_to_seed_ratio;
		}
	}
	map<string, ProbModel> trait_distributions;
	map<int, string> distribution_types = { { 0, "uniform" }, {1, "linear"}, {2, "normal"}, {3, "discrete"}, {4, "constant"} };
};


class Tree {
public:
	Tree() = default;
	Tree(int _id, pair<float, float> _position, float _dbh, float seed_bearing_threshold,
		map<int, float> _resprout_growthcurve, float _growth_multiplier
	) :
		position(_position), dbh(_dbh), resprout_growthcurve(_resprout_growthcurve)
	{
		id = _id;
		derive_allometries(seed_bearing_threshold);
		growth_multiplier = _growth_multiplier;
		//printf("Tree created with id: %i, radius: %f, radius_tmin1: %f, stem dbh: %f, bark thickness: %f, LAI: %f \n", id, radius, radius_tmin1, dbh, bark_thickness, LAI);
	};
	bool operator==(const Tree& tree) const
	{
		return id == tree.id;
	}
	void resprout(float seed_bearing_threshold) {
		life_phase = 1;
		dbh = 0;
		age = 0;
		derive_allometries(seed_bearing_threshold);
	}
	bool derive_allometries(float seed_bearing_threshold) {
		crown_area = compute_crown_area();
		radius = compute_radius();
		auto [_life_phase, became_reproductive] = get_life_phase(seed_bearing_threshold);
		life_phase = _life_phase;
		bark_thickness = get_bark_thickness();
		LAI = get_LAI();
		height = get_height();
		lowest_branch = get_lowest_branch_height();
		return became_reproductive;
	}
	bool radius_spans(pair<float, float> pos2, bool verbose = false) {
		float dist = help::get_dist(position, pos2);
		if (verbose) printf("	Distance between tree %i and given position %f, %f: %f \n", id, pos2.first, pos2.second, dist);
		return dist < radius;
	}
	pair<int, bool> get_life_phase(float& seed_bearing_threshold) {
		int new_life_phase;
		bool life_phase_changed = false;
		if (dbh > seed_bearing_threshold) new_life_phase = 2;
		else new_life_phase = life_phase;
		if (new_life_phase != life_phase) {
			life_phase_changed = true;
		}
		return pair<int, bool>(new_life_phase, life_phase_changed);
	}
	float get_bark_thickness() {
		return 0.31 * pow(dbh, 1.276); // From Hoffmann et al (2012), figure 5a. Bark thickness in mm.
	}
	float get_dbh_from_radius() {
		float basal_area = pow(10.0f, (log10(crown_area) + 0.32f) / 0.59f); // From Rossatto et al (2009), figure 6. Reordered equation.
		return sqrtf(basal_area / M_PI); // Convert basal area (cm^2) to dbh (cm).
	}
	float get_dbh_increment(float LAI_shade) {
		if (LAI_shade > 5.0f) return 0.0f; // If the LAI of shading leaf cover is larger than 5, the tree is too shaded to grow.

		float stem_increment = 0.3f * (1.0f - exp(-0.118 * dbh * 10.0f));	// I = I_max(1-e^(-g * D)), from Hoffman et al (2012), Appendix 1, page 1.
																	// Current stem dbh and increment in cm.

		stem_increment *= (5.0f - LAI_shade) / 5.0f;	// LAI-dependent growth reduction to introduce density dependence. Multiplication factor: ((LAI_max - LAI) / LAI_max) 
														// From Hoffman et al (2012), Appendix 2.


		return stem_increment;
	}
	float get_survival_probability(float& fire_resistance_argmin, float& fire_resistance_argmax, float& fire_resistance_stretch) {
		return help::get_sigmoid(bark_thickness, fire_resistance_argmin, fire_resistance_argmax, fire_resistance_stretch); // Sigmoid based on Hoffman et al (2012), figure 2a.
	}
	float get_leaf_area() {
		return 0.147 * pow(dbh, 2.053); // From Hoffman et al (2012), figure 5b. Leaf area in m^2
	}
	float compute_crown_area() {
		return 0.551f * pow(dbh, 1.28f);	// Crown area in m^2. y = b x^a. Parameters of a and b are averages of fits for 5 trees presented in Blanchard et al (2015),
											// page 1957 (explains the model in "Fitting tree allometries"), and table 4 (shows regression results).
	}
	float get_LAI() {
		float leaf_area = get_leaf_area();
		return leaf_area / crown_area;
	}
	bool survives_fire(float &fire_resistance_argmin, float &fire_resistance_argmax, float &fire_resistance_stretch) {
		float survival_probability = get_survival_probability(fire_resistance_argmin, fire_resistance_argmax, fire_resistance_stretch);
		return help::get_rand_float(0.0f, 1.0f) < survival_probability;
	}
	float compute_radius() {
		float radius_breast_height = dbh * 0.5f; // Convert dbh (cm) to stem radius at breast height (cm).
		float basal_area = M_PI * radius_breast_height * radius_breast_height;
		crown_area = pow(10.0f, 0.59 * log10(basal_area) - 0.32); // From Rossatto et al (2009), figure 6.
		float crown_radius = sqrt(crown_area / M_PI);
		return crown_radius;
	}
	float get_height() {
		float ln_dbh = log(dbh);
		return exp(0.865 + 0.760 * ln_dbh - 0.0340 * (ln_dbh * ln_dbh)); // From Chave et al (2014), equation 6a. Value of E obtained by calculating mean of E
																		 // values for bistable study sites.
	}
	float get_lowest_branch_height() {
		return height * 0.4f; // We assume the tree's crown begins at 40% its height.
							  // TODO: Perhaps make this fraction a function of dbh for added realism.
	}
	float get_AGB() {
		float ln_dbh = log(dbh);
		float ln_wood_specific_gravity = log(0.5f);
		return -1.803 - 0.976f * -0.02802 + 0.976 * ln_wood_specific_gravity + 2.673f * ln_dbh - 0.0299f * (ln_dbh * ln_dbh); // From Chave et al (2014), equation 7.
	}
	float compute_new_dbh(float LAI_shade) {
		float _dbh;
		if (dbh < 2.5f) {
			if (life_phase == 1) {
				_dbh = resprout_growthcurve.at(age); // Resprouts younger than 5 years (implied by dbh < 2.5) are assumed to grow according to a predefined growth curve (Hoffmann et al, 2012, supplementary information 1).
			}
			else {
				_dbh = dbh + growth_multiplier * 0.25f; // Assume a constant growth rate of 2.5 mm for saplings, until they reach 2.5 cm dbh. Based on growth rate of resprouts after first 5 years (Hoffmann et al, 2012, supplementary information 1. Also see "Tree Allometric Relations.xlsx").
			}
		}
		else {
			_dbh = dbh + growth_multiplier * get_dbh_increment(LAI_shade);
		}
		return _dbh;
	}
	pair<bool, bool> grow(float &seed_bearing_threshold, float LAI_shade) {
		age++;
		float _dbh = compute_new_dbh(LAI_shade);
		bool dies_due_to_light_limitation = is_float_equal(_dbh, dbh) && (life_phase == 0); // If the tree is not reproductive yet and is unable to grow, we assume it dies.
		dbh = _dbh;
		bool became_reproductive = derive_allometries(seed_bearing_threshold);
		return pair<bool, bool>(became_reproductive, dies_due_to_light_limitation);
	}
	void print() {
		printf("id: %d, dbh: %f, radius: %f, bark thickness: %f, LAI: %f, height: %f, lowest branch: %f, crown area: %f, position: %f, %f \n",
			id, dbh, radius, bark_thickness, LAI, height, lowest_branch, crown_area, position.first, position.second
		);
	}
	float radius = -1;
	float dbh = 0.1;
	float bark_thickness = 0;
	float shade = 0;
	float LAI = 0;
	float height = 0;
	float lowest_branch = 0;
	float crown_area = 0;
	float growth_multiplier = 1;
	pair<float, float> position = pair(0, 0);
	int id = -1;
	int age = -1;
	int life_phase = 0;
	int last_mortality_check = 0;
	map<int, float> resprout_growthcurve;
};


class Crop {
public:
	Crop() = default;
	Crop(Strategy &_strategy, Tree& tree) {
		strategy = _strategy;
		seed_mass = _strategy.seed_mass;
		origin = tree.position;
		id = tree.id;
	}
	void compute_no_seeds(Tree &tree, float STR) {
		float dbh_dependent_factor = (tree.dbh / 30.0f);
		total_no_seeds_produced = STR * (dbh_dependent_factor * dbh_dependent_factor); // Number of seeds produced by the tree, based on Ribbens et al (1994).
		no_seeds = (float)total_no_seeds_produced * strategy.recruitment_probability;
	}
	void compute_no_diaspora() {
		no_diaspora = no_seeds / strategy.no_seeds_per_diaspore;
		no_seeds = no_diaspora * strategy.no_seeds_per_diaspore; // Ensure that the number of seeds is an exact multiple of the number of seeds per diaspore.
	}
	void compute_fruit_abundance() {
		// Compute the number of fruits that will be produced by the tree, rather than the number of fruits that are actually dispersed.
		fruit_abundance = total_no_seeds_produced / strategy.no_seeds_per_diaspore;
	}
	void update(Tree& tree, float STR) {
		compute_no_seeds(tree, STR);
		compute_no_diaspora();
		compute_fruit_abundance();
	}
	int no_seeds = 0; // Number of seeds dispersed.
	int no_diaspora = 0; // Number of diaspora dispersed.
	int fruit_abundance = 0;
	int total_no_seeds_produced = 0; // Number of seeds produced (always equal- or higher than no_seeds).
	float seed_mass = 0;
	Strategy strategy;
	pair<float, float> origin = pair<float, float>(0, 0);
	int id = -1;
};


class Population {
public:
	Population() = default;
	Population(float _max_dbh, float _cellsize, float dbh_q1, float dbh_q2,
		map<string, map<string, float>> strategy_parameters, float _mutation_rate, float _seed_bearing_threshold,
		float growth_multiplier_stdev, float growth_multiplier_min, float growth_multiplier_max
	) : max_dbh(_max_dbh), cellsize(_cellsize), seed_bearing_threshold(_seed_bearing_threshold)
	{
		strategy_generator = StrategyGenerator(strategy_parameters);
		dbh_probability_model = help::LinearProbabilityModel(dbh_q1, dbh_q2, 0, max_dbh);
		mutation_rate = _mutation_rate;
		init_growth_curves();
		growth_multiplier_distribution = ProbModel(1, growth_multiplier_stdev, growth_multiplier_min, growth_multiplier_max, -1);
	}
	void init_growth_curves() {
		resprout_growthcurve = {
			{1, 1.3f}, {2, 1.8f}, {3, 2.1f}, {4, 2.5f}, // From Hoffmann et al (2012), estimated from supplementary figure S1.
		};
	}
	Tree* add(pair<float, float> position, Strategy* _strategy = 0, float dbh = -2) {
		// Create tree
		if (dbh == -1) dbh = max_dbh;
		else if (dbh == -2) {
			if (_strategy == nullptr) {
				while (dbh <= 0)
					dbh = dbh_probability_model.linear_sample();
			}
			else {
				dbh = _strategy->seedling_dbh; // Growth rate determines initial dbh.
			}
		}
		float growth_multiplier = help::get_rand_float(growth_multiplier_distribution.min_value, growth_multiplier_distribution.max_value);
		Tree tree(no_created_trees + 1, position, dbh, seed_bearing_threshold, resprout_growthcurve, growth_multiplier);
		members[tree.id] = tree;
		no_created_trees++;

		// Create strategy
		Strategy strategy;
		if (_strategy != nullptr) {
			strategy = *_strategy;
			strategy_generator.mutate(strategy, mutation_rate);
		}
		else {
			strategy_generator.generate(strategy);
		} 
		strategy.id = tree.id;
		recruitment_rates.push_back(strategy.recruitment_probability);

		// Create crop
		Crop crop(strategy, tree);
		crops[tree.id] = crop;

		// Create custom kernel
		Kernel kernel = *get_kernel(strategy.vector);
		if (strategy.vector == "wind") {
			kernel = Kernel(
				tree.id, kernel.dist_max, kernel.wspeed_gmean, kernel.wspeed_stdev, kernel.wind_direction,
				kernel.wind_direction_stdev, strategy.seed_tspeed
			);
		}
		kernels_individual[tree.id] = kernel;

		return &members[tree.id];
	}
	void add_reproduction_system(Tree &tree) {

	}
	Kernel* add_kernel(string tree_dispersal_vector, Kernel &kernel) {
		kernels[tree_dispersal_vector] = kernel;
		return &kernel;
	}
	Tree* get(int id) {
		return &members[id];
	}
	void get(vector<int>& ids, vector<Tree*> &trees) {
		for (int id : ids) trees.push_back(get(id));
	}
	Crop* get_crop(int id) {
		return &crops[id];
	}
	Kernel* get_kernel(string vector) {
		return &kernels[vector];
	}
	Kernel* get_kernel(int id) {
		return &kernels_individual[id];
	}
	int size() {
		return members.size();
	}
	bool remove(Tree* tree) {
		return remove(tree->id);
	}
	bool remove(int id) {
		bool removed = members.erase(id);
		removed = removed && crops.erase(id);
		removed = removed && delete_kernel(id);
		return removed;
	}
	bool delete_kernel(int id) {
		if (get_kernel(id)->type == "wind") delete[] kernels_individual[id].cdf;
		return kernels_individual.erase(id);
	}
	bool is_population_member(Tree* tree) {
		auto it = members.find(tree->id);
		return (it != members.end());
	}
	bool is_population_member(int tree_id) {
		return members.find(tree_id) != members.end();
	}
	void free() {
		for (auto& [id, tree] : members) {
			delete_kernel(id);
		}
	}
	void get_ids_and_trait_values(string trait, map<int, double> &ids_and_trait_values) {
		if (trait == "height") {
			for (auto& [id, tree] : members) {
				ids_and_trait_values[id] = tree.height;
			}
		}
	}
	void sort_by_trait(string trait, PairSet &sorted_population) {
		map<int, double> ids_and_trait_values;
		get_ids_and_trait_values(trait, ids_and_trait_values);
		help::sort(ids_and_trait_values, sorted_population);
	}
	
	unordered_map<int, Tree> members;
	unordered_map<int, Crop> crops;
	unordered_map<string, Kernel> kernels;
	unordered_map<int, Kernel> kernels_individual;
	help::LinearProbabilityModel dbh_probability_model;
	StrategyGenerator strategy_generator;
	float max_dbh = 0;
	float cellsize = 0;
	float seed_mass = 0;
	float mutation_rate = 0;
	vector<float> recruitment_rates;
	Tree removed_tree;
	int no_created_trees = 0;
	float seed_bearing_threshold = 0;
	map<int, float> resprout_growthcurve;
	ProbModel growth_multiplier_distribution;
};
