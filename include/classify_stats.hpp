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

double variance(const auto & v, const auto & mean) {
    auto variance_func = [&mean](float accumulator, const float& val) {
        return accumulator + ((val - mean)*(val - mean));
    };
    double sum = std::accumulate(std::begin(v), std::end(v), 0.0, variance_func);
    return v.size() <= 1 ? 0 : sum / (v.size() - 1);
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
            PLOG_VERBOSE << "Add to g_pos results in g_pos size " << pos.size() << " and g_neg size " << neg.size() << " with status " << pos_complete << neg_complete << complete;
            return complete;
        };

        bool add_neg(const float& val){
            if (neg.size() < num_reads_to_fit and val > 0){
#pragma omp critical(neg_push)
                neg.push_back(val);
            } else {
                check_status();
            }
            PLOG_VERBOSE << "Add to g_neg results in g_pos size " << pos.size() << " and g_neg size " << neg.size() << " with status " << pos_complete << neg_complete << complete;
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

    void fit_loc(const std::vector<float> & training_data)
    {
        const auto mu = mean(training_data);
        loc = mu - (shape * scale);
    }

    double calculate_anderson_darling(std::vector<float> & training_data)
    {
        // in ascending order
        sort(training_data.begin(), training_data.end());

        auto n = training_data.size();
        double S = 0;
        for (auto i=0; i<training_data.size(); ++i)
        {
            auto F_y_i = stats::pgamma(training_data[i]-loc,shape, scale);
            auto F_y_n1i = stats::pgamma(training_data[n-1-i]-loc,shape, scale);
            auto s_i = ((2*i + 1)/n) * (std::log(F_y_i) + std::log(1-F_y_n1i));
            S += s_i;
        }
        return -n-S;
    }
};

struct BetaParams
{
    float alpha;
    float beta;
    float loc;

    BetaParams(const float & _alpha, const float & _beta):
            alpha{_alpha},
            beta{_beta}
    {};

    void fit_loc(const std::vector<float> & training_data)
    {
        const auto mu = mean(training_data);
        loc = mu - alpha/(alpha + beta);
    }

    void fit (const std::vector<float> & training_data, bool is_pos)
    {
        const auto mu = mean(training_data);
        const auto var = variance(training_data, mu);
        PLOG_INFO << "Fitting Beta to data with mean " << mu << " and sample variance " << var;
        assert(var < mu*(1-mu));

        alpha = mu*((mu*(1-mu)/var)-1);
        beta = (1-mu)*((mu*(1-mu)/var)-1);

        // rough way to force support of negative beta to cover the region [0,0.2] - larger numbers narrow the supported peak
        if (beta > 85)
            beta = 85;

        // rough way to force pos distribution to be right skewed
        if (is_pos and beta > alpha)
            alpha = beta;

        assert(alpha > 0);
        assert(beta > 0);
    }

    double calculate_anderson_darling(std::vector<float> & training_data)
    {
        // in ascending order
        sort(training_data.begin(), training_data.end());

        auto n = training_data.size();
        double S = 0;
        for (auto i=0; i<training_data.size(); ++i)
        {
            auto F_y_i = stats::pbeta(training_data[i],alpha, beta);
            auto F_y_n1i = stats::pbeta(training_data[n-1-i],alpha, beta);
            auto s_i = ((2*i + 1)/(n-1)) * (std::log(F_y_i) + std::log(1-F_y_n1i));
            S += s_i;
        }
        return -n+1-S;
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
        std::string dist{"gamma"};
        GammaParams g_pos{25, 0, 0.02};
        GammaParams g_neg{10, 0, 0.005};
        BetaParams b_pos{6,4};
        BetaParams b_neg{5,80};
    public:
        Model() = default;
        Model(Model const &) = default;
        Model(Model &&) = default;
        Model & operator=(Model const &) = default;
        Model & operator=(Model &&) = default;
        ~Model() = default;

        Model(const uint8_t i, const std::string& d):
            id{i},
            dist{d}
        {};

        void train_gamma(TrainingData & training_data)
        {
            assert(training_data.id == id);
            if (training_data.pos_complete ) {
                g_pos.fit(training_data.pos);
                auto ad = g_pos.calculate_anderson_darling(training_data.pos);
                PLOG_INFO << "Model " << +id << " fit g_pos data with Gamma (shape:" << g_pos.shape << ", loc: 0, scale: "
                          << g_pos.scale << "). Anderson-darling statistic is " << ad;
            } else {
                auto ad = g_pos.calculate_anderson_darling(training_data.pos);
                PLOG_INFO << "Model " << +id << " using default for g_pos data with Gamma (shape:" << g_pos.shape
                          << ", loc: 0, scale: " << g_pos.scale << "). Anderson-darling statistic is " << ad;
            }

            if (training_data.neg_complete ) {
                g_neg.fit(training_data.neg);
                auto ad = g_neg.calculate_anderson_darling(training_data.neg);
                PLOG_INFO << "Model " << +id << " fit g_neg data with Gamma (shape:" << g_neg.shape << ", loc: 0, scale: "
                          << g_neg.scale << "). Anderson-darling statistic is " << ad;
            } else {
                g_neg.fit_loc(training_data.neg);
                auto ad = g_neg.calculate_anderson_darling(training_data.neg);
                PLOG_INFO << "Model " << +id << " using default for g_neg data with Gamma (shape:" << g_neg.shape
                          << ", loc: " << g_neg.loc << ", scale: " << g_neg.scale << "). Anderson-darling statistic is " << ad;
            }
        }

        void train_beta(TrainingData & training_data)
        {
            assert(training_data.id == id);
            if (training_data.pos_complete ) {
                b_pos.fit(training_data.pos, true);
                auto ad = b_pos.calculate_anderson_darling(training_data.pos);
                PLOG_INFO << "Model " << +id << " fit pos data with Beta (alpha:" << b_pos.alpha << ", beta: "
                          << b_pos.beta << "). Anderson-darling statistic is " << ad;
            } else {
                auto ad = b_pos.calculate_anderson_darling(training_data.pos);
                PLOG_INFO << "Model " << +id << " using default for pos data with Beta (alpha:" << b_pos.alpha << ", beta: "
                          << b_pos.beta << "). Anderson-darling statistic is " << ad;
            }

            if (training_data.neg_complete ) {
                b_neg.fit(training_data.neg, false);
                auto ad = b_neg.calculate_anderson_darling(training_data.neg);
                PLOG_INFO << "Model " << +id << " fit neg data with Beta (alpha:" << b_neg.alpha << ", beta: "
                          << b_neg.beta << ", loc: " << b_neg.loc << "). Anderson-darling statistic is " << ad;
            } else {
                auto ad = b_neg.calculate_anderson_darling(training_data.neg);
                PLOG_INFO << "Model " << +id << " using default for neg data with Beta (alpha:" << b_neg.alpha << ", beta: "
                          << b_neg.beta << ", loc: " << b_neg.loc << "). Anderson-darling statistic is " << ad;
            }
        }

        void train(TrainingData & training_data)
        {
            assert(training_data.id == id);
            if (dist == "beta")
                train_beta(training_data);
            else
                train_gamma(training_data);
            ready = true;
            training_data.clear();
        }

        ProbPair prob(const float & read_proportion) const {
            const auto p_err = stats::dexp(read_proportion,300);
            float p_pos, p_neg;
            if ( dist == "gamma") {
                p_pos = stats::dgamma(read_proportion - g_pos.loc, g_pos.shape, g_pos.scale);
                p_neg = stats::dgamma(read_proportion - g_neg.loc, g_neg.shape, g_neg.scale);
            } else {
                p_pos = stats::dbeta(read_proportion, b_pos.alpha, b_pos.beta);
                p_neg = stats::dbeta(read_proportion - b_neg.loc, b_neg.alpha, b_neg.beta);
            }
            if (read_proportion == 1)
                p_pos = 1;

            const auto total = p_err + p_pos + p_neg;
            // Use Neyman Pearson Lemma
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
        uint8_t min_hits_;

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
                confidence_threshold_(opt.confidence_threshold),
                min_hits_(opt.min_hits)
        {
            for (auto i=0; i<summary.num_categories(); ++i){
                models_.emplace_back(Model(i, opt.dist));
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

        uint8_t min_num_hits() const
        {
            return min_hits_;
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
            PLOG_VERBOSE << "add_to_training is " << add_to_training << " with hi g_pos " << +pos_i;

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
