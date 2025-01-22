#ifndef __CLASSIFY_STATS_H_INCLUDED__
#define __CLASSIFY_STATS_H_INCLUDED__

#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cassert>

#include <stats.hpp>

#include "input_summary.hpp"
#include "classify_arguments.hpp"

double mean(const auto & v) {
    double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
    return v.empty() ? 0 : sum / v.size();
}

class TrainingData
{
    private:
        uint8_t id;
        bool complete = false;
        bool pos_complete = false;
        bool neg_complete = false;
        uint16_t num_reads_to_fit;
        std::vector<float> pos;
        std::vector<float> neg;
    public:
        TrainingData() = default;
        TrainingData(TrainingData const &) = default;
        TrainingData(TrainingData &&) = default;
        TrainingData & operator=(TrainingData const &) = default;
        TrainingData & operator=(TrainingData &&) = default;
        ~TrainingData() = default;

        TrainingData(const ClassifyArguments & opt, const uint8_t i):
            id{i}
        {
            num_reads_to_fit = opt.num_reads_to_fit;
            pos.reserve(opt.num_reads_to_fit);
            neg.reserve(opt.num_reads_to_fit);
        };

        bool check_status()
        {
            if (pos.size() >= num_reads_to_fit)
                pos_complete = true;
            if (neg.size() >= num_reads_to_fit)
                neg_complete = true;
            if (pos_complete and neg_complete)
                complete = true;
            return complete;
        }

        bool add_pos(const float& val){
            if (pos.size() < num_reads_to_fit){
#pragma omp critical(pos_push)
                pos.push_back(val);
            } else {
                check_status();
            }
            PLOG_VERBOSE << "Add to pos results in pos size " << pos.size() << " and neg size " << neg.size() << " with status " << pos_complete << neg_complete << complete;
            return complete;
        };

        bool add_neg(const float& val){
            if (neg.size() < num_reads_to_fit and val > 0){
#pragma omp critical(neg_push)
                neg.push_back(val);
            } else {
                check_status();
            }
            PLOG_VERBOSE << "Add to neg results in pos size " << pos.size() << " and neg size " << neg.size() << " with status " << pos_complete << neg_complete << complete;
            return complete;
        };

        void clear()
        {
            pos.clear();
            pos.resize(0);
            neg.clear();
            neg.resize(0);
        }

    friend class StatsModel;
    friend class Model;
};

struct GammaParams
{
    float shape;
    float loc;
    float scale;

    GammaParams(const float & _shape, const float & _loc, const float & _scale):
        shape{_shape},
        loc{_loc},
        scale{_scale}
    {};

    void fit (const std::vector<float> & training_data)
    {
        const auto mu = mean(training_data);
        const auto ln_mu = std::log(mu);
        const auto log_data = training_data | std::views::transform([](const auto& n) { return std::log(n); });

        const auto mean_ln = mean(log_data);
        const auto s = ln_mu - mean_ln;
        shape = (3 - s + std::sqrt((s-3)*(s-3) + 24*s))/(12*s);
        scale = mu / shape;
    }
};

struct ProbPair {
    double pos;
    double neg;
};

class Model
{
    private:
        uint8_t id;
        bool ready {false};
        GammaParams pos{25,0,0.02};
        GammaParams neg{10,0,0.005};
    public:
        Model() = default;
        Model(Model const &) = default;
        Model(Model &&) = default;
        Model & operator=(Model const &) = default;
        Model & operator=(Model &&) = default;
        ~Model() = default;

        Model(const uint8_t i):
            id{i}
        {};

        void train(TrainingData & training_data)
        {
            assert(training_data.id == id);
            if (training_data.pos_complete ) {
                pos.fit(training_data.pos);
                PLOG_INFO << "Model " << +id << " fit pos data with Gamma (shape:" << pos.shape <<", loc: 0, scale: "
                          << pos.scale << ")";
            } else {
                PLOG_INFO << "Model " << +id << " using default for pos data with Gamma (shape:" << pos.shape <<", loc: 0, scale: "<< pos.scale << ")";
            }

            if (training_data.neg_complete ) {
                neg.fit(training_data.neg);
                PLOG_INFO << "Model " << +id << " fit neg data with Gamma (shape:" << neg.shape << ", loc: 0, scale: "
                          << neg.scale << ")";
            } else {
                PLOG_INFO << "Model " << +id << " using default for neg data with Gamma (shape:" << neg.shape <<", loc: 0, scale: "<< neg.scale << ")";
            }
            ready = true;
            training_data.clear();
        }

        ProbPair prob(const float & read_proportion) const {
            const auto p_err = stats::dexp(read_proportion,1000);
            const auto p_pos = stats::dgamma(read_proportion,pos.shape,pos.scale);
            const auto p_neg = stats::dgamma(read_proportion,neg.shape,neg.scale);
            const auto total = p_err + p_pos + p_neg;
            return ProbPair(p_pos/total, (p_err + p_neg)/total);
        }

    friend class StatsModel;

};

class StatsModel
{
    private:
        bool ready_ {false};
        float lo_hi_threshold_;
        int8_t confidence_threshold_;
        std::vector<TrainingData> training_data_;
        std::vector<Model> models_;

    public:
        StatsModel() = default;
        StatsModel(StatsModel const &) = default;
        StatsModel(StatsModel &&) = default;
        StatsModel & operator=(StatsModel const &) = default;
        StatsModel & operator=(StatsModel &&) = default;
        ~StatsModel() = default;

        StatsModel(const ClassifyArguments& opt, const InputSummary & summary):
                lo_hi_threshold_(opt.lo_hi_threshold),
                confidence_threshold_(opt.confidence_threshold)
        {
            for (auto i=0; i<summary.num_categories(); ++i){
                models_.emplace_back(Model(i));
                training_data_.emplace_back(TrainingData(opt, i));
            }
            PLOG_DEBUG << "Initialize stats model with " << training_data_.size() << " sets of training data ";
        };

        bool ready() const
        {
            return ready_;
        }

        void force_ready()
        {
            for (auto i=0; i < models_.size(); ++i){
                auto & model = models_[i];
                if (model.ready)
                    continue;
                auto & data = training_data_[i];
#pragma omp critical(train)
                model.train(data);
            }
            ready_ = true;
        }

        void check_if_ready()
        {
            if (ready_)
                return;
            for (const auto & model : models_){
                if (not model.ready)
                    return;
            }
            ready_ = true;
        }

        int8_t confidence_threshold() const
        {
            return confidence_threshold_;
        }

    void train_model_at(const uint8_t & i)
        {
            PLOG_DEBUG << "Train model at position " << +i;
            auto & data = training_data_[i];
            assert(data.complete and data.pos_complete and data.neg_complete);

            auto & model = models_[i];
            assert(not model.ready);
#pragma omp critical(train)
            {
                model.train(data);
                check_if_ready();
            }
        }

        bool add_read_to_training_data(const std::vector<float>& read_proportions){
            auto pos_i = std::numeric_limits<uint8_t>::max();
            auto max_val = 0.0;
            auto num_above_threshold = 0;
            PLOG_VERBOSE << "Decide if read is training candidate";
            for (uint8_t i=0; i < read_proportions.size(); ++i) {
                const auto & val = read_proportions.at(i);
                PLOG_VERBOSE << "Read prop " << val << " at " << +i;
                if (val > lo_hi_threshold_)
                    num_above_threshold += 1;
                if (val == max_val)
                    pos_i = std::numeric_limits<uint8_t>::max();
                else if (val > max_val)
                {
                    pos_i = i;
                    max_val = val;
                }
            }
            bool add_to_training = (pos_i != std::numeric_limits<uint8_t>::max() and num_above_threshold <= 1);
            PLOG_VERBOSE << "add_to_training is " << add_to_training << " with hi pos " << +pos_i;

            if (add_to_training) {
                PLOG_VERBOSE << " read_proportions size is " << read_proportions.size() << " and training data partition has size " << training_data_.size();
                auto ready_to_train = training_data_.at(pos_i).add_pos(read_proportions[pos_i]);
                if (ready_to_train and not models_[pos_i].ready){
#pragma omp critical(train_model)
                    train_model_at(pos_i);
                }

                for (uint8_t i=0; i < read_proportions.size(); ++i) {
                    if (i != pos_i){
                        ready_to_train = training_data_.at(i).add_neg(read_proportions[i]);
                        if (ready_to_train and not models_[i].ready) {
#pragma omp critical(train_model)
                            train_model_at(i);
                        }
                    }
                }
            }
            return ready_;
        };

        ProbPair classify(const auto & i, const float & read_proportion) const
        {
            const auto & model = models_.at(i);
            return model.prob(read_proportion);
        }
};

#endif
