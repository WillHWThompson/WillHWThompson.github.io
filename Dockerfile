FROM ruby:3.1
LABEL MAINTAINER Will Thompson

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libxml2-dev \
    imagemagick \
    git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /srv/jekyll

# Copy Gemfile and install dependencies
COPY Gemfile* /srv/jekyll/
# Install gems into vendor/bundle and use bundler v2
ENV BUNDLE_PATH=/srv/jekyll/vendor/bundle
RUN bundle install --jobs 4 --retry 3

# Copy site files into working directory
ADD . /srv/jekyll/

# Set up environment
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Set bundler to install gems into vendor/bundle
ENV BUNDLE_PATH=/srv/jekyll/vendor/bundle
