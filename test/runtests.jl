using TestItemRunner

if get(ENV, "GITHUB_ACTIONS", "false") == "true"
    @run_package_tests verbose=true filter=ti->!(:skip_ci in ti.tags)
else
    @run_package_tests verbose=true
end
