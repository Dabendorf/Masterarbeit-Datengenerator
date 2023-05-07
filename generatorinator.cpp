#include <cmath>
#include <cassert>
#include <cstring>
#include <sys/time.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <vector>
#include <algorithm>
#include <typeinfo>

#define FEATUREVECTORS_FILE "gen/chr22_feature.vectors"
#define GENES_FILE "gen/genes.txt"
#define SKYSERVER_DATA "sky/skyserver.txt"
#define SKYSERVER_PREDEFINED_QUERIES "sky/skyserver.queries"

static double gettime(void) {
	struct timeval now_tv;
	gettimeofday (&now_tv,NULL);
	return ((double)now_tv.tv_sec) + ((double)now_tv.tv_usec) / 1000000.0;
}

int main(int argc, char* argv[]) {
	size_t  n;
	size_t  m;
	double   o;
	double sel;
	int rq;
	int experiment_type;
	int query_type;

	if (argc < 5) {
		std::cout << "Requires at least 5 arguments" << std::endl;
		return 1;

	} else {
		n = atoi(argv[1]);
		m = atoi(argv[2]);
		o = 1.0;
		sel = 0.2;
		experiment_type = atoi(argv[3]);
		rq = atoi(argv[4]);
	}
	std::vector< std::vector<double> > data_points;

	if(experiment_type==7) {
		data_points = std::vector< std::vector<double> >(n, std::vector<double>(m+1));
	} else {
		data_points = std::vector< std::vector<double> >(n, std::vector<double>(m));
	}

	bool is_gmrqb;
	bool is_skyserver;

	if (experiment_type == 1 || experiment_type == 2 || experiment_type == 4 || experiment_type == 8 || experiment_type == 9) {
		is_gmrqb = true;
		is_skyserver = false;
	} else {
		is_gmrqb = false;
		is_skyserver = true;
	}

	if (is_gmrqb) {
		size_t i = 0;
		std::ifstream feature_vectors(FEATUREVECTORS_FILE);
		std::string line;
		std::string token;

		while (std::getline(feature_vectors, line) && i < n) {
			std::vector<double> data_point(m);
			std::vector<std::string> line_tokens;
			std::istringstream iss(line);
			while(std::getline(iss, token, ' ')) {
				line_tokens.push_back(token);
			}

			for (size_t j = 0; j < m; ++j) {
				data_point[j] = stod(line_tokens[j]);
			}
			data_points[i++] = data_point;
		}
	} else if(is_skyserver && experiment_type != 7) {
		size_t i = 0;
		std::ifstream feature_vectors(SKYSERVER_DATA);
		std::string line;
		std::string token;
		std::getline(feature_vectors, line);

		while (std::getline(feature_vectors, line) && i < n) {
			std::vector<double> data_point(m);
			std::vector<std::string> line_tokens;
			std::istringstream iss(line);
			while(std::getline(iss, token, ',')) {
				line_tokens.push_back(token);
			}
				
			for (size_t j = 0; j < m; ++j) {
				data_point[j] = stod(line_tokens[j]);
			}
			data_points[i++] = data_point;
		}
	} else if(is_skyserver && experiment_type == 7) {
		size_t i = 0;
		std::ifstream feature_vectors(SKYSERVER_DATA);
		std::string line;
		std::string token;
		std::getline(feature_vectors, line);
		std::vector< std::vector<double> > data_points_temp(n*(m/2)+m%2, std::vector<double>(3));

		while (std::getline(feature_vectors, line) && i < n*(m/2)+m%2) {
			std::vector<double> data_point(3);
			std::vector<std::string> line_tokens;
			std::istringstream iss(line);
			while(std::getline(iss, token, ',')) {
				line_tokens.push_back(token);
			}
				
			for (size_t j = 0; j < 3; ++j) {
				data_point[j] = stod(line_tokens[j]);
			}
			data_points_temp[i++] = data_point;
		}

		auto rng = std::default_random_engine {};
		std::shuffle(std::begin(data_points_temp), std::end(data_points_temp), rng);

		for(int ind_datapoint_temp=0; ind_datapoint_temp < n; ind_datapoint_temp++) {
			
			data_points[ind_datapoint_temp][0] = ind_datapoint_temp;
			int temp_count = 1;
			for(int s=0; s<(m/2)+m%2; s++) {
				for(int i=1; i<3; i++) {
					data_points[ind_datapoint_temp][temp_count++] = data_points_temp[ind_datapoint_temp*(m/2)+m%2+s][i];
				}
			}
		}
	}

	std::ofstream data_outputfile;
	if (is_gmrqb) {
		data_outputfile.open("gmrqb_datapoints.txt");
	} else {
		data_outputfile.open("skyserver_datapoints.txt");
	}

	for (unsigned int i = 0; i < data_points.size(); ++i) {
		for (unsigned int j = 0; j < data_points[i].size(); ++j) {
				data_outputfile << data_points[i][j] << " ";
		}
		data_outputfile << std::endl;
	}

	data_outputfile.close();
	
	std::ofstream query_outputfile;
	if (is_gmrqb) {
		query_outputfile.open("gmrqb_range.queries");
	} else {
		query_outputfile.open("skyserver_range.queries");
	}
	
	std::vector<std::vector<double> > lb_queries(rq, std::vector<double>(m, std::numeric_limits<double>::min()));
	std::vector<std::vector<double> > ub_queries(rq, std::vector<double>(m, std::numeric_limits<double>::max()));

	if (experiment_type == 1 || experiment_type == 2 || experiment_type == 9) {
		size_t i = 0;
		std::ifstream genes(GENES_FILE);
		std::string line;
		std::string token;

		while (std::getline(genes, line) && i < rq) {
			std::vector<std::string> line_tokens;
			std::istringstream iss(line);
			while(std::getline(iss, token, '\t')) {
				line_tokens.push_back(token);
			}
				
			lb_queries[i][5] = (double) 22;
			ub_queries[i][5] = (double) 22;
			lb_queries[i][6] = (double) std::stod(line_tokens[4]) - 100000.0;
			ub_queries[i][6] = (double) std::stod(line_tokens[5]) + 200000.0;
			
			if (experiment_type == 1) {
				if (argc >= 6) {
					query_type = atoi(argv[5]);
				} else {
					std::cout << "missing argument" << std::endl;
					return 1;
				}
			} else {
				query_type = rand() % 8;
			}

			int rand_point = rand() % n;

			switch (query_type) {
				// Query 2
				case 1:
					// qual (create range around a certain qual found in the data set)
					lb_queries[i][8] = data_points[rand_point][8] * 0.5;
					ub_queries[i][8] = lb_queries[i][8] * 3;
					// depth (create range around a certain depth found in the data set)
					lb_queries[i][9] = data_points[rand_point][9] * 0.5;
					ub_queries[i][9] = lb_queries[i][9] * 3;
					// allele freq (create range using a certain allele_freq found in the data set)
					lb_queries[i][10] = data_points[rand_point][10];
					ub_queries[i][10] = lb_queries[i][10] + 0.3;
					break;
				// Query 3
				case 2:
					// gender
					lb_queries[i][2] = data_points[rand_point][2];
					ub_queries[i][2] = lb_queries[i][2];
					break;
				// Query 4 
				case 3:
					// gender
					lb_queries[i][2] = data_points[rand_point][2];
					ub_queries[i][2] = lb_queries[i][2];
					// population
					lb_queries[i][1] = data_points[rand_point][1];
					ub_queries[i][1] = lb_queries[i][1];
					break;
				// Query 5
				case 4:
					// gender
					lb_queries[i][2] = data_points[rand_point][2];
					ub_queries[i][2] = lb_queries[i][2];
					// population
					lb_queries[i][1] = data_points[rand_point][1];
					ub_queries[i][1] = lb_queries[i][1];
					// relationship
					lb_queries[i][4] = data_points[rand_point][4];
					ub_queries[i][4] = lb_queries[i][4];
					break;
				// Query 6
				case 5:
					// gender
					lb_queries[i][2] = data_points[rand_point][2];
					ub_queries[i][2] = lb_queries[i][2];
					// population
					lb_queries[i][1] = data_points[rand_point][1];
					ub_queries[i][1] = lb_queries[i][1];
					// relationship
					lb_queries[i][4] = data_points[rand_point][4];
					ub_queries[i][4] = lb_queries[i][4];
					// family_id (create range using a certain family_id found in the data set)
					lb_queries[i][3] = data_points[rand_point][3] * 0.5;
					ub_queries[i][3] = lb_queries[i][3] * 3;
					break;
				// Query 7
				case 6:
					// gender
					lb_queries[i][2] = data_points[rand_point][2];
					ub_queries[i][2] = lb_queries[i][2];
					// population
					lb_queries[i][1] = data_points[rand_point][1];
					ub_queries[i][1] = lb_queries[i][1];
					// relationship
					lb_queries[i][4] = data_points[rand_point][4];
					ub_queries[i][4] = lb_queries[i][4];
					// family_id (create range using a certain family_id found in the data set)
					lb_queries[i][3] = data_points[rand_point][3] * 0.5;
					ub_queries[i][3] = lb_queries[i][3] * 3;
					// mutation_id (create range using a certain mutation_id found in the data set)
					lb_queries[i][7] = data_points[rand_point][7] * 0.5;
					ub_queries[i][7] = lb_queries[i][7] * 3;
					break;
				// Query 8
				case 7:
					// sample_id
					lb_queries[i][0] = data_points[rand_point][0];
					ub_queries[i][0] = lb_queries[i][0];
					// population
					lb_queries[i][1] = data_points[rand_point][1];
					ub_queries[i][1] = lb_queries[i][1];
					// gender
					lb_queries[i][2] = data_points[rand_point][2];
					ub_queries[i][2] = lb_queries[i][2];
					// family_id (create range using a certain family_id found in the data set)
					lb_queries[i][3] = data_points[rand_point][3] * 0.5;
					ub_queries[i][3] = lb_queries[i][3] * 3;
					// relationship
					lb_queries[i][4] = data_points[rand_point][4];
					ub_queries[i][4] = lb_queries[i][4];
					// mutation_id (create range using a certain mutation_id found in the data set)
					lb_queries[i][7] = data_points[rand_point][7] * 0.5;
					ub_queries[i][7] = lb_queries[i][7] * 3;
					// qual (create range around a certain qual found in the data set)
					lb_queries[i][8] = data_points[rand_point][8] * 0.5;
					ub_queries[i][8] = lb_queries[i][8] * 3;
					// depth (create range around a certain depth found in the data set)
					lb_queries[i][9] = data_points[rand_point][9] * 0.5;
					ub_queries[i][9] = lb_queries[i][9] * 3;
					// allele freq (create range using a certain allele_freq found in the data set)
					lb_queries[i][10] = data_points[rand_point][10];
					ub_queries[i][10] = lb_queries[i][10] + 0.3;
					// allele_count (create range using a certain allele_count found in the data set)
					lb_queries[i][11] = data_points[rand_point][11] * 0.5;
					ub_queries[i][11] = lb_queries[i][11] * 3;
					// filter
					lb_queries[i][12] = data_points[rand_point][12];
					ub_queries[i][12] = lb_queries[i][12];
					// ref_base
					lb_queries[i][13] = data_points[rand_point][13];
					ub_queries[i][13] = lb_queries[i][13];
					// alt_base
					lb_queries[i][14] = data_points[rand_point][14];
					ub_queries[i][14] = lb_queries[i][14];
					// variant_type
					lb_queries[i][15] = data_points[rand_point][15];
					ub_queries[i][15] = lb_queries[i][15];
					// ancestral_allele
					lb_queries[i][16] = data_points[rand_point][16];
					ub_queries[i][16] = lb_queries[i][16];
					// genotypegenotype
					lb_queries[i][17] = data_points[rand_point][17];
					ub_queries[i][17] = lb_queries[i][17];
					// reference genome
					lb_queries[i][18] = data_points[rand_point][18];
					ub_queries[i][18] = lb_queries[i][18];
					break;
				// Query 1
				case 0:
				default:
					break;
			}
			i++;
		}
	} else if(experiment_type == 3) {
		size_t i = 0;
		std::ifstream genes(SKYSERVER_PREDEFINED_QUERIES);
		std::string line;
		std::string token;

		while (std::getline(genes, line) && i < rq) {
			std::vector<std::string> line_tokens;
			std::istringstream iss(line);
			while(std::getline(iss, token, ',')) {
				line_tokens.push_back(token);
			}

			query_type = rand() % 3;
			std::cout << "Query type: " << query_type << std::endl;

			int rand_point = rand() % n;

			switch (query_type) {
				case 0:
					lb_queries[i][0] = std::stod(line_tokens[0]);
					ub_queries[i][0] = std::stod(line_tokens[1]);
					break;
				case 1:
					lb_queries[i][1] = std::stod(line_tokens[2]);
					ub_queries[i][1] = std::stod(line_tokens[3]);
					break;
				case 2:
					lb_queries[i][0] = std::stod(line_tokens[0]);
					ub_queries[i][0] = std::stod(line_tokens[1]);
					lb_queries[i][1] = std::stod(line_tokens[2]);
					ub_queries[i][1] = std::stod(line_tokens[3]);
					break;
				default:
					break;
			}
			i++;
		}
	} else if(experiment_type == 7) {
		size_t i = 0;
		std::ifstream genes(SKYSERVER_PREDEFINED_QUERIES);
		std::string line;
		std::string token;

		while (std::getline(genes, line) && i < rq) {
			std::vector<std::string> line_tokens;
			for(int s=0; s<(m/2)+m%2; s++) {
				std::istringstream iss(line);
				while(std::getline(iss, token, ',')) {
					line_tokens.push_back(token);
				}
			}

			for(int s=0; s<m; s++) {
				lb_queries[i][s] = std::stod(line_tokens[2*s]);
				ub_queries[i][s] = std::stod(line_tokens[2*s+1]);
			}
			i++;
		}
	} else if(experiment_type == 4) {
		// used similar code to cracking kd project
		lb_queries = std::vector<std::vector<double>>();
		ub_queries = std::vector<std::vector<double>>();

		int workload_type;
		if (argc >= 6) {
			workload_type = atoi(argv[5]);
		} else {
			std::cout << "missing argument" << std::endl;
			return 1;
		}

		std::vector<std::vector<double>> sorted_points(n, std::vector<double>(m));

		for (size_t col = 0; col < m; ++col) {
			std::vector<double> column(n);
			for (size_t row = 0; row < n; ++row) {
				column[row] = data_points[row][col];
			}

			std::sort(column.begin(), column.end());

			for (size_t row = 0; row < n; ++row) {
				sorted_points[row][col] = column[row];
			}
		}

		size_t n_rows = n;
		size_t n_dimensions = m;

		if(workload_type==1) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			std::mt19937 generator_query(1);
			std::uniform_int_distribution<int> distr_query(
				0, n_rows*(1-per_column_sel)
			);

			for(size_t i = 0; i < rq; ++i){
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = distr_query(generator_query);
					highs.at(j) = lows.at(j) + (int64_t)(n_rows * per_column_sel);
					cols.at(j) = j;
				}
				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		} else if(workload_type==2) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			double mean = n_rows/2.0;
			double std_dev = mean/8.0;

			std::mt19937 generator_query(1);
			std::normal_distribution<double> distr_query(mean, std_dev);

			for(size_t i = 0; i < rq; ++i){
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					auto v = distr_query(generator_query);
					lows.at(j) = v - (per_column_sel/2.0 * n_rows);
					highs.at(j) = v + (per_column_sel/2.0 * n_rows);
					cols.at(j) = j;
				}
				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		} else if(workload_type==3) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			double side = (n_rows * per_column_sel)/5.0;
			double half_side = side/2.0;

			auto center = half_side;

			for(size_t i = 0; i < rq && center <= n_rows - side; ++i, center += side) {
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = center - half_side;
					highs.at(j) = center + half_side;
					cols.at(j) = j;
				}
				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		} else if(workload_type==5) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			auto center = n_rows/2.0;

			double half_side = (n_rows * per_column_sel)/2.0;
			
			double step = half_side / rq;

			for(size_t i = rq; i > 0; --i){
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = center - i*step;
					highs.at(j) = center + i*step;
					cols.at(j) = j;
				}

				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		} else if(workload_type==6) {
			sel = 1.0;
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);
			double half_side = (n_rows * per_column_sel)/2.0;
			auto number_of_blocks = n_rows/(2*half_side);
			size_t rq_per_block = rq/number_of_blocks;
			double step = half_side / rq_per_block;

			for(size_t i = 0; i < number_of_blocks; ++i){
				double center = 2 * half_side * (i + 1);
				for(size_t s = 0 ; s < rq_per_block; ++s){
					std::vector<double> lows(n_dimensions);
					std::vector<double> highs(n_dimensions);
					std::vector<size_t> cols(n_dimensions);

					for(size_t j = 0; j < n_dimensions; ++j){
						lows.at(j) = (center-half_side) + s*step;
						highs.at(j) = (center+half_side) - s*step;
						cols.at(j) = j;
					}
					std::vector<double> lb_temp(19);
					std::vector<double> ub_temp(19);
					for (unsigned int j = 0; j < lows.size(); ++j) {
						lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
						ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
					}
					lb_queries.push_back(lb_temp);
					ub_queries.push_back(ub_temp);
					}
			}
		} else if(workload_type==7) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			double half_side = (n_rows * per_column_sel)/2.0;

			float center_low = half_side;
			float center_high = n_rows - half_side;
			
			double step = half_side / rq;

			for(size_t i = rq; i > rq/2; --i){
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = center_low - i*step;
            		highs.at(j) = center_low + i*step;
					cols.at(j) = j;
				}
				
				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = center_high - i*step;
					highs.at(j) = center_high + i*step;
					cols.at(j) = j;
				}

				lb_temp.assign(m, 0.0);
				ub_temp.assign(m, 0.0);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		}

	} else if(experiment_type == 5) {
		// used similar code to cracking kd project
		lb_queries = std::vector<std::vector<double>>();
		ub_queries = std::vector<std::vector<double>>();

		int workload_type;
		if (argc >= 6) {
			workload_type = atoi(argv[5]);
		} else {
			std::cout << "missing argument" << std::endl;
			return 1;
		}

		std::vector<std::vector<double>> sorted_points(n, std::vector<double>(m));

		for (size_t col = 0; col < m; ++col) {
			std::vector<double> column(n);
			for (size_t row = 0; row < n; ++row) {
				column[row] = data_points[row][col];
			}

			std::sort(column.begin(), column.end());

			for (size_t row = 0; row < n; ++row) {
				sorted_points[row][col] = column[row];
			}
		}

		size_t n_rows = n;
		size_t n_dimensions = m;

		if(workload_type==1) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			std::mt19937 generator_query(1);
			std::uniform_int_distribution<int> distr_query(
				0, n_rows*(1-per_column_sel)
			);

			for(size_t i = 0; i < rq; ++i){
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = distr_query(generator_query);
					highs.at(j) = lows.at(j) + (int64_t)(n_rows * per_column_sel);
					cols.at(j) = j;
				}
				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		} else if(workload_type==2) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			double mean = n_rows/2.0;
			double std_dev = mean/8.0;

			std::mt19937 generator_query(1);
			std::normal_distribution<double> distr_query(mean, std_dev);

			for(size_t i = 0; i < rq; ++i){
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					auto v = distr_query(generator_query);
					lows.at(j) = v - (per_column_sel/2.0 * n_rows);
					highs.at(j) = v + (per_column_sel/2.0 * n_rows);
					cols.at(j) = j;
				}
				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		} else if(workload_type==3) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			double side = (n_rows * per_column_sel)/5.0;
			double half_side = side/2.0;

			auto center = half_side;

			for(size_t i = 0; i < rq && center <= n_rows - side; ++i, center += side) {
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = center - half_side;
					highs.at(j) = center + half_side;
					cols.at(j) = j;
				}
				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		} else if(workload_type==5) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			auto center = n_rows/2.0;

			double half_side = (n_rows * per_column_sel)/2.0;
			
			double step = half_side / rq;

			for(size_t i = rq; i > 0; --i){
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = center - i*step;
					highs.at(j) = center + i*step;
					cols.at(j) = j;
				}

				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		} else if(workload_type==6) {
			sel = 1.0;
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);
			double half_side = (n_rows * per_column_sel)/2.0;
			auto number_of_blocks = n_rows/(2*half_side);
			size_t rq_per_block = rq/number_of_blocks;
			double step = half_side / rq_per_block;

			for(size_t i = 0; i < number_of_blocks; ++i){
				double center = 2 * half_side * (i + 1);
				for(size_t s = 0 ; s < rq_per_block; ++s){
					std::vector<double> lows(n_dimensions);
					std::vector<double> highs(n_dimensions);
					std::vector<size_t> cols(n_dimensions);

					for(size_t j = 0; j < n_dimensions; ++j){
						lows.at(j) = (center-half_side) + s*step;
						highs.at(j) = (center+half_side) - s*step;
						cols.at(j) = j;
					}
					std::vector<double> lb_temp(19);
					std::vector<double> ub_temp(19);
					for (unsigned int j = 0; j < lows.size(); ++j) {
						lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
						ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
					}
					lb_queries.push_back(lb_temp);
					ub_queries.push_back(ub_temp);
					}
			}
		} else if(workload_type==7) {
			double per_column_sel = std::pow(sel, 1.0/n_dimensions);

			double half_side = (n_rows * per_column_sel)/2.0;

			float center_low = half_side;
			float center_high = n_rows - half_side;
			
			double step = half_side / rq;

			for(size_t i = rq; i > rq/2; --i){
				std::vector<double> lows(n_dimensions);
				std::vector<double> highs(n_dimensions);
				std::vector<size_t> cols(n_dimensions);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = center_low - i*step;
            		highs.at(j) = center_low + i*step;
					cols.at(j) = j;
				}
				
				std::vector<double> lb_temp(19);
				std::vector<double> ub_temp(19);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);

				for(size_t j = 0; j < n_dimensions; ++j){
					lows.at(j) = center_high - i*step;
					highs.at(j) = center_high + i*step;
					cols.at(j) = j;
				}

				lb_temp.assign(m, 0.0);
				ub_temp.assign(m, 0.0);
				for (unsigned int j = 0; j < lows.size(); ++j) {
					lb_temp[j] = sorted_points[std::max(lows[j],(double)0)][j];
					ub_temp[j] = sorted_points[std::min(highs[j],(double)(n_rows-1))][j];
				}
				lb_queries.push_back(lb_temp);
				ub_queries.push_back(ub_temp);
			}
		}

	} else if(experiment_type == 6) {
		double selectivity;
		if (argc >= 6) {
			selectivity = std::stod(argv[5]);
		} else {
			std::cout << "missing argument" << std::endl;
			return 1;
		}

		double sel_per_dim;
		if(m==3) {
			sel_per_dim = sqrt(selectivity*0.01);
		} else {
			sel_per_dim = std::pow(selectivity, 1.0/(m-1));
		}
		
		if(m==3) {
			std::vector<double> data_points_dim1;
			std::vector<double> data_points_dim2;
			data_points_dim1 = std::vector<double>(n);
			data_points_dim2 = std::vector<double>(n);

			for(int i=0; i < n; i++) {
				data_points_dim1[i] = data_points[i][1];
				data_points_dim2[i] = data_points[i][2];
			}

			std::sort(data_points_dim1.begin(), data_points_dim1.end());
			std::sort(data_points_dim2.begin(), data_points_dim2.end());
			
			for(int i=0; i<rq; i++) {
				int lower_bound = rand() % (int) (n*(1-sel_per_dim));
				int upper_bound = std::min((int) (n-1), (int) (lower_bound+n*sel_per_dim));
				lb_queries[i][0] = data_points_dim1[lower_bound];
				ub_queries[i][0] = data_points_dim1[upper_bound];

				lower_bound = rand() % (int)(n*(1-sel_per_dim));
				upper_bound = std::min((int) (n-1), (int) (lower_bound+n*sel_per_dim));
				lb_queries[i][1] = data_points_dim2[lower_bound];
				ub_queries[i][1] = data_points_dim2[upper_bound];
			}
		}
	}

	if(is_gmrqb) {
		if(experiment_type != 8 && experiment_type != 4) {
			for (size_t i = 0; i < rq; ++i) {
				if (lb_queries[i][0] == std::numeric_limits<double>::min()) {
					query_outputfile << "min";
				} else {
					query_outputfile << lb_queries[i][0];
				}

				for (size_t j = 1; j < m; ++j) {
					if (lb_queries[i][j] == std::numeric_limits<double>::min()) {
						query_outputfile << " min";
					}else{
						query_outputfile << " " << lb_queries[i][j];
					}
				}
				query_outputfile << std::endl;

				if (ub_queries[i][0] == std::numeric_limits<double>::max()) {
					query_outputfile << "max";
				}
				else {
					query_outputfile << ub_queries[i][0];
				}

				for (size_t j = 1; j < m; ++j) {
					if (ub_queries[i][j] == std::numeric_limits<double>::max()) {
						query_outputfile << " max";
					} else {
						query_outputfile << " " << ub_queries[i][j];
					}
				}
				query_outputfile << std::endl;
			}
		} else if(experiment_type == 4) {
			for (size_t i = 0; i < rq; ++i) {
				if (lb_queries[i][0] == std::numeric_limits<double>::min()) {
					query_outputfile << "min";
				} else {
					query_outputfile << lb_queries[i][0];
				}

				for (size_t j = 1; j < 19; ++j) {
					if (lb_queries[i][j] == std::numeric_limits<double>::min() || lb_queries[i][j] == 0) {
						query_outputfile << " min";
					}else{
						query_outputfile << " " << lb_queries[i][j];
					}
				}
				query_outputfile << "\n";

				if (ub_queries[i][0] == std::numeric_limits<double>::max()) {
					query_outputfile << "max";
				}
				else {
					query_outputfile << ub_queries[i][0];
				}

				for (size_t j = 1; j < 19; ++j) {
					if (ub_queries[i][j] == std::numeric_limits<double>::max() || ub_queries[i][j] == 0) {
						query_outputfile << " max";
					} else {
						query_outputfile << " " << ub_queries[i][j];
					}
				}
				query_outputfile << "\n";
			}
		}
	} else if(is_skyserver) {
		if(experiment_type != 8) {
			int num_dim;
			if(experiment_type==7) {
				num_dim = m;
			} else {
				num_dim = 2;
			}
			for (size_t i = 0; i < rq; ++i) {
				if (lb_queries[i][0] == std::numeric_limits<double>::min()) {
					query_outputfile << "min";
				} else {
					query_outputfile << lb_queries[i][0];
				}

				for (size_t j = 1; j < num_dim; ++j) {
					if (lb_queries[i][j] == std::numeric_limits<double>::min()) {
						query_outputfile << " min";
					}else{
						query_outputfile << " " << lb_queries[i][j];
					}
				}
				query_outputfile << std::endl;

				if (ub_queries[i][0] == std::numeric_limits<double>::max()) {
					query_outputfile << "max";
				}
				else {
					query_outputfile << ub_queries[i][0];
				}

				for (size_t j = 1; j < num_dim; ++j) {
					if (ub_queries[i][j] == std::numeric_limits<double>::max()) {
						query_outputfile << " max";
					} else {
						query_outputfile << " " << ub_queries[i][j];
					}
				}
				query_outputfile << std::endl;
			}
		}
	}

	query_outputfile.close();

	if(experiment_type == 8 || experiment_type == 9) {
		std::vector<std::vector<double>> points_to_query(0, std::vector<double>(m));

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dist(0, data_points.size() - 1);
		std::vector<int> indices(data_points.size());
		std::iota(indices.begin(), indices.end(), 0);
		std::shuffle(indices.begin(), indices.end(), gen);
		for (int i = 0; i < rq; ++i) {
			points_to_query.push_back(data_points[indices[i]]);
		}

		std::ofstream pointsquery_outputfile;
		pointsquery_outputfile.open("gmrqb_points.queries");

		for (unsigned int i = 0; i < points_to_query.size(); ++i) {
			for (unsigned int j = 0; j < points_to_query[i].size(); ++j) {
					pointsquery_outputfile << points_to_query[i][j] << " ";
			}
			pointsquery_outputfile << std::endl;
		}

		pointsquery_outputfile.close();
	}

	return 0;
}
