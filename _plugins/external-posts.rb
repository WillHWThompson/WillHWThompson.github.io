require 'feedjira'
require 'httparty'
require 'jekyll'

module ExternalPosts
  class ExternalPostsGenerator < Jekyll::Generator
    safe true
    priority :high

    def generate(site)
      if site.config['external_sources'] != nil
        site.config['external_sources'].each do |src|
          begin
            puts "Fetching external posts from #{src['name']}:"
            response = HTTParty.get(src['rss_url'], timeout: 10)
            
            if response.code != 200
              puts "Error fetching RSS feed from #{src['rss_url']} - HTTP Status: #{response.code}"
              next
            end
            
            begin
              xml = response.body
              feed = Feedjira.parse(xml)
              
              feed.entries.each do |e|
                puts "...fetching #{e.url}"
                slug = e.title.downcase.strip.gsub(' ', '-').gsub(/[^\w-]/, '')
                path = site.in_source_dir("_posts/#{slug}.md")
                doc = Jekyll::Document.new(
                  path, { :site => site, :collection => site.collections['posts'] }
                )
                doc.data['external_source'] = src['name'];
                doc.data['feed_content'] = e.content;
                doc.data['title'] = "#{e.title}";
                doc.data['description'] = e.summary;
                doc.data['date'] = e.published;
                doc.data['redirect'] = e.url;
                site.collections['posts'].docs << doc
              end
            rescue Feedjira::NoParserAvailable => e
              puts "Error parsing feed from #{src['rss_url']}: #{e.message}"
              puts "Feed content (first 500 chars):"
              puts xml.to_s[0..500]
            rescue => e
              puts "Unexpected error processing feed from #{src['rss_url']}: #{e.message}"
              puts e.backtrace.join("\n")
            end
          rescue HTTParty::Error, SocketError, Timeout::Error, Errno::ECONNREFUSED => e
            puts "Network error fetching feed from #{src['rss_url']}: #{e.message}"
          rescue => e
            puts "Unexpected error while fetching from #{src['rss_url']}: #{e.class} - #{e.message}"
          end
        end
      end
    end
  end
end
